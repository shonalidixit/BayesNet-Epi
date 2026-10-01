"""BayesNet-Epi: reproducible synthetic network epidemiology demonstration."""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'figures'; DATA=ROOT/'data'

def sigmoid(x): return 1/(1+np.exp(-x))

def simulate_network(n=120, seed=42):
    rng=np.random.default_rng(seed)
    regions=rng.integers(0,4,n)
    G=nx.DiGraph(); G.add_nodes_from(range(n))
    for i in range(n):
        G.nodes[i]['region']=int(regions[i])
    for i in range(n):
        for j in range(n):
            if i==j: continue
            p=.065 if regions[i]==regions[j] else .012
            if rng.random()<p:
                G.add_edge(i,j,weight=int(rng.poisson(4)+1))
    return G

def simulate_outbreak(G, steps=10, beta=.055, seed=42):
    rng=np.random.default_rng(seed); nodes=list(G.nodes)
    infected={int(rng.choice(nodes))}; infection_time={i:(0 if i in infected else -1) for i in nodes}
    for t in range(1,steps+1):
        new=set()
        for u in list(infected):
            for _,v,d in G.out_edges(u,data=True):
                if v not in infected and rng.random()<1-math.exp(-beta*d['weight']): new.add(v)
        for v in new: infection_time[v]=t
        infected |= new
    return infection_time

def features(G, infection_time):
    pr=nx.pagerank(G,weight='weight')
    rows=[]
    for n in G.nodes:
        rows.append({'farm':n,'region':G.nodes[n]['region'],'in_degree':G.in_degree(n),'out_degree':G.out_degree(n),
                     'weighted_in':G.in_degree(n,weight='weight'),'weighted_out':G.out_degree(n,weight='weight'),
                     'pagerank':pr[n],'infected':int(infection_time[n]>=0),'infection_time':infection_time[n]})
    return pd.DataFrame(rows)

def fit_bayesian_logistic_laplace(df, seed=42):
    """Bayesian logistic regression via MAP + Laplace approximation, no specialist PPL required."""
    from scipy.optimize import minimize
    X=df[['weighted_in','weighted_out','pagerank']].to_numpy(float)
    X=(X-X.mean(0))/(X.std(0)+1e-9); X=np.c_[np.ones(len(X)),X]; y=df.infected.to_numpy(float)
    prior_sd=2.5
    def nlp(w):
        z=X@w
        ll=np.sum(y*z-np.logaddexp(0,z))
        return -ll + .5*np.sum((w/prior_sd)**2)
    res=minimize(nlp,np.zeros(X.shape[1]),method='BFGS')
    w=res.x; p=sigmoid(X@w); W=p*(1-p)
    H=X.T@(X*W[:,None])+np.eye(X.shape[1])/(prior_sd**2)
    cov=np.linalg.inv(H); se=np.sqrt(np.diag(cov))
    names=['Intercept','Weighted in-degree','Weighted out-degree','PageRank']
    return pd.DataFrame({'parameter':names,'posterior_mode':w,'approx_sd':se,'lower_95':w-1.96*se,'upper_95':w+1.96*se})

def community_analysis(G):
    U=G.to_undirected()
    comm=list(nx.algorithms.community.greedy_modularity_communities(U,weight='weight'))
    mapping={n:i for i,c in enumerate(comm) for n in c}
    true=np.array([G.nodes[n]['region'] for n in G.nodes]); pred=np.array([mapping[n] for n in G.nodes])
    # adjusted Rand without sklearn
    from scipy.special import comb
    table=pd.crosstab(true,pred).to_numpy(); n=len(true)
    sum_comb=sum(comb(x,2) for x in table.ravel()); a=sum(comb(x,2) for x in table.sum(1)); b=sum(comb(x,2) for x in table.sum(0)); total=comb(n,2)
    expected=a*b/total; ari=(sum_comb-expected)/(.5*(a+b)-expected) if .5*(a+b)!=expected else 0
    return mapping,float(ari)

def sampling_experiment(G, fractions=(.25,.5,.75), reps=20, seed=42):
    rng=np.random.default_rng(seed); base=nx.transitivity(G.to_undirected()); rows=[]; nodes=np.array(G.nodes)
    for f in fractions:
        for r in range(reps):
            chosen=rng.choice(nodes,size=max(3,int(len(nodes)*f)),replace=False); H=G.subgraph(chosen).copy()
            rows.append({'fraction':f,'rep':r,'density':nx.density(H),'transitivity':nx.transitivity(H.to_undirected()),'transitivity_error':abs(nx.transitivity(H.to_undirected())-base)})
    return pd.DataFrame(rows)

def plot_network(G, infection_time):
    pos=nx.spring_layout(G,seed=7); vals=[infection_time[n] if infection_time[n]>=0 else 12 for n in G.nodes]
    plt.figure(figsize=(9,7)); nx.draw_networkx_edges(G,pos,alpha=.12,arrows=False); nx.draw_networkx_nodes(G,pos,node_size=35,node_color=vals,cmap='viridis'); plt.axis('off'); plt.title('Synthetic directed farm-movement network\n(node shade reflects infection timing; uninfected assigned final shade)'); plt.tight_layout(); plt.savefig(FIG/'network_outbreak.png',dpi=180); plt.close()

def main():
    FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)
    G=simulate_network(); it=simulate_outbreak(G); df=features(G,it); post=fit_bayesian_logistic_laplace(df); mapping,ari=community_analysis(G); samp=sampling_experiment(G)
    df['community']=df.farm.map(mapping); df.to_csv(DATA/'synthetic_farms.csv',index=False); post.to_csv(DATA/'bayesian_posterior.csv',index=False); samp.to_csv(DATA/'sampling_results.csv',index=False)
    nx.write_gexf(G,DATA/'movement_network.gexf'); plot_network(G,it)
    summary={'farms':G.number_of_nodes(),'movements':G.number_of_edges(),'infected':int(df.infected.sum()),'community_ARI_vs_simulated_regions':ari,'mean_sampling_transitivity_error':float(samp.transitivity_error.mean())}
    (DATA/'summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2)); print('\nPosterior:\n',post.to_string(index=False))
if __name__=='__main__': main()

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pipeline import simulate_network, simulate_outbreak, features

def test_network_is_directed_and_nonempty():
    G=simulate_network(30,1); assert G.is_directed(); assert G.number_of_edges()>0

def test_outbreak_and_features_cover_nodes():
    G=simulate_network(30,2); t=simulate_outbreak(G,seed=2); df=features(G,t); assert len(t)==30; assert len(df)==30; assert set(df.infected.unique()) <= {0,1}

import unittest
import os
import pandas as pd
import networkx as nx 
from src.main import parse_pdb, process_dataframe, plot_atom_point_cloud, create_graphein_object, plot_backbone_graph

class TestProteinStructureVisualizer(unittest.TestCase):

    def setUp(self):
        self.test_pdb_path = "test.pdb"
        self.mock_pdb_content = """
ATOM      1  N   MET A   1      20.154  34.103  23.570  1.00  0.00           N  
ATOM      2  CA A MET A   1      21.268  33.204  23.870  1.00  0.00           C  
ATOM      3  CA B MET A   1      21.568  33.504  23.870  1.00  0.00           C  
ATOM      4  C   MET A   1      21.961  33.613  25.223  1.00  0.00           C  
ATOM      5  O   MET A   1      21.479  34.572  25.857  1.00  0.00           O  
ATOM      6  CB A MET A   1      22.316  33.285  22.764  1.00  0.00           C  
ATOM      7  CB B MET A   1      22.616  33.585  22.764  1.00  0.00           C  
TER
END
        """
        with open(self.test_pdb_path, "w") as f:
            f.write(self.mock_pdb_content)
        
        self.output_path = "test.output.html"
    
    def tearDown(self):
        if os.path.exists(self.test_pdb_path):
            os.remove(self.test_pdb_path)
        if os.path.exists(self.output_path):
            os.remove(self.output_path)
    
    def test_parse_pdb(self):
        atomic_df, header = parse_pdb(self.test_pdb_path)
        self.assertIsInstance(atomic_df, pd.DataFrame)
        self.assertGreater(len(atomic_df), 0)
    
    def test_process_dataframe(self):
        atomic_df, _ = parse_pdb(self.test_pdb_path)
        processed_df = process_dataframe(atomic_df)
        self.assertIsInstance(processed_df, pd.DataFrame)
        self.assertTrue("node_id" in processed_df.columns)
        self.assertTrue(processed_df["atom_name"].eq("CA").all())
    
        self.assertTrue(processed_df["alt_loc"].isin(["A", ""]).all())
        self.assertFalse(processed_df["alt_loc"].eq("B").any())
    
    def test_plot_atom_point_cloud(self):
        atomic_df, _ = parse_pdb(self.test_pdb_path)
        plot_atom_point_cloud(atomic_df, self.output_path)
        self.assertTrue(os.path.exists(self.output_path))
    
    def test_create_graphein_object(self):
        atomic_df, _ = parse_pdb(self.test_pdb_path)
        processed_df = process_dataframe(atomic_df)
        graphein_object = create_graphein_object(processed_df, atomic_df, "test")
        self.assertIsInstance(graphein_object, nx.Graph)
        self.assertGreater(len(graphein_object.nodes), 0)
    
    def test_plot_backbone_graph(self):
        atomic_df, _ = parse_pdb(self.test_pdb_path)
        processed_df = process_dataframe(atomic_df)
        graphein_object = create_graphein_object(processed_df, atomic_df, "test")
        plot_backbone_graph(graphein_object, "test", self.output_path)
        self.assertTrue(os.path.exists(self.output_path))

if __name__ == "__main__":
    unittest.main()

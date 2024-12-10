import unittest
import os
import pandas as pd
import networkx as nx
from src.main import (
    parse_pdb,
    remove_altloc,
    create_graphein_object,
    plot_atom_point_cloud,
    plot_backbone_graph,
    DEFAULT_GRANULARITY,
    OUTPUT_FILES
)

class TestProteinStructureVisualizer(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create test directories
        self.output_dir = "test_output"
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create test PDB file
        self.test_pdb_path = os.path.join(self.output_dir, "test.pdb")
        self.mock_pdb_content = """
ATOM      1  N   MET A   1      20.154  34.103  23.570  1.00 30.00           N  
ATOM      2  CA  MET A   1      21.268  33.204  23.870  1.00 30.00           C  
ATOM      3  C   MET A   1      21.961  33.613  25.223  1.00 30.00           C  
ATOM      4  CA  MET A   2      22.268  34.204  24.870  1.00 30.00           C  
ATOM      5  C   MET A   2      22.961  34.613  26.223  1.00 30.00           C  
END
"""
        with open(self.test_pdb_path, "w") as f:
            f.write(self.mock_pdb_content)
        
        # Set up output paths
        self.point_cloud_path = os.path.join(self.output_dir, "test_point_cloud.html")
        self.backbone_path = os.path.join(self.output_dir, "test_backbone.html")
        
        self.sample_df = pd.DataFrame({
            'atom_name': ['CA', 'CA', 'CA'],
            'alt_loc': ['', '', ''],
            'residue_name': ['MET', 'MET', 'MET'],
            'chain_id': ['A', 'A', 'A'],
            'residue_number': [1, 2, 3],
            'x_coord': [1.0, 2.0, 3.0],
            'y_coord': [1.0, 2.0, 3.0],
            'z_coord': [1.0, 2.0, 3.0],
            'occupancy': [1.0, 1.0, 1.0],
            'b_factor': [30.0, 30.0, 30.0],
            'element_symbol': ['C', 'C', 'C'],
            'charge': ['', '', ''],
            'line_idx': [0, 1, 2],
            'model_id': [1, 1, 1]
        })

    def tearDown(self):
        """Clean up test fixtures after each test method."""
        if os.path.exists(self.test_pdb_path):
            os.remove(self.test_pdb_path)
        
        for file in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, file))
        os.rmdir(self.output_dir)

    def test_parse_pdb(self):
        """Test PDB file parsing."""
        df, header = parse_pdb(self.test_pdb_path)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)
        self.assertTrue(all(col in df.columns for col in ['atom_name', 'x_coord', 'y_coord', 'z_coord']))

    def test_parse_pdb_invalid_file(self):
        """Test parsing an invalid PDB file."""
        with self.assertRaises(FileNotFoundError):
            parse_pdb("nonexistent.pdb")

    def test_remove_altloc(self):
        """Test alternative location removal."""
        processed_df = remove_altloc(self.sample_df)
        self.assertIsInstance(processed_df, pd.DataFrame)
        self.assertTrue(all(processed_df['atom_name'] == DEFAULT_GRANULARITY))

    def test_remove_altloc_invalid_input(self):
        """Test remove_altloc with invalid input."""
        # Test with string input
        with self.assertRaises(TypeError) as context:
            remove_altloc("not_a_dataframe")
        self.assertEqual(str(context.exception), "Input must be a pandas DataFrame")
        
        # Test with None
        with self.assertRaises(TypeError) as context:
            remove_altloc(None)
        self.assertEqual(str(context.exception), "Input must be a pandas DataFrame")

    def test_create_graphein_object(self):
        """Test Graphein object creation."""
        processed_df = remove_altloc(self.sample_df)
        graph = create_graphein_object(processed_df, self.sample_df, "test.pdb")
        self.assertIsInstance(graph, nx.Graph)
        self.assertGreater(len(graph.nodes), 0)

    def test_plot_atom_point_cloud(self):
        """Test point cloud plot generation."""
        plot_atom_point_cloud(self.sample_df, self.point_cloud_path)
        self.assertTrue(os.path.exists(self.point_cloud_path))
        self.assertGreater(os.path.getsize(self.point_cloud_path), 0)

    def test_plot_backbone_graph(self):
        """Test backbone graph plot generation."""
        processed_df = remove_altloc(self.sample_df)
        graph = create_graphein_object(processed_df, self.sample_df, "test.pdb")
        plot_backbone_graph(graph, "test", self.backbone_path)
        self.assertTrue(os.path.exists(self.backbone_path))
        self.assertGreater(os.path.getsize(self.backbone_path), 0)

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        empty_df = pd.DataFrame(columns=self.sample_df.columns)
        with self.assertRaises(ValueError) as context:
            remove_altloc(empty_df)
        self.assertEqual(str(context.exception), "Input DataFrame is empty")

    def test_missing_columns(self):
        """Test handling of DataFrame with missing required columns."""
        incomplete_df = self.sample_df.drop(columns=['alt_loc'])
        with self.assertRaises(KeyError) as context:
            remove_altloc(incomplete_df)
        self.assertIn("Missing required columns", str(context.exception))

if __name__ == '__main__':
    unittest.main()

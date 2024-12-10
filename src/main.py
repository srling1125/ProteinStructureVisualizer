import pandas as pd
import plotly.express as px
import networkx as nx
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import webbrowser
import os
from biopandas.pdb import PandasPdb
from prody import parsePDBHeader
from typing import Optional, Tuple
from graphein.protein.visualisation import plotly_protein_structure_graph
from graphein.protein.graphs import label_node_id
from graphein.protein.graphs import initialise_graph_with_metadata
from graphein.protein.graphs import add_nodes_to_graph

# Constants and Configurations
DEFAULT_MODEL_INDEX = 1
DEFAULT_GRANULARITY = "CA"
OUTPUT_FILES = {
    "point_cloud": "point_cloud.html",
    "backbone_graph": "backbone_graph.html"
}

def parse_pdb(pdb_path: Optional[str] = None, model_index: int = 1, parse_header: bool = True) -> pd.DataFrame:
    """
    Parses a PDB file into a dataframe using PandasPdb. 

    Parameters:
    - pdb_path (str): Path to the PDB file.
    - model_index (int): Model index to extract (default: 1).
    - parse_header (bool): Whether to parse the PDB header (default: True). 

    Returns: 
    - pd.DataFrame: Concatenated dataframe of ATOM and HETATM data. 
    - dict: Parsed PDB header or None if parse_header is False.
    """
    parsed_df = PandasPdb().read_pdb(pdb_path)
    if parse_header:
        header = parsePDBHeader(pdb_path)
    else:
        header = None
    parsed_df = parsed_df.get_model(model_index)
    if len(parsed_df.df["ATOM"]) == 0:
        raise ValueError(f"No model found for index: {model_index}")
    
    return pd.concat([parsed_df.df["ATOM"], parsed_df.df["HETATM"]]), header

def remove_altloc(df: pd.DataFrame, granularity: str = DEFAULT_GRANULARITY) -> pd.DataFrame:
    """
    Processes the parsed dataframe to remove alternative locations of atoms.
    
    Parameters: 
    - df (pd.DataFrame): Input atomic data (parsed_df) from parse_pdb. 
    - granularity (str): Granularity for graph (default is "CA" for alpha carbon)
    
    Returns:
    - pd.DataFrame: Processed DataFrame.
    
    Raises:
    - TypeError: If input is not a pandas DataFrame
    - ValueError: If DataFrame is empty
    - KeyError: If required columns are missing
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame")
        
    if df.empty:
        raise ValueError("Input DataFrame is empty")
        
    required_columns = ['atom_name', 'alt_loc', 'chain_id']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise KeyError(f"Missing required columns: {missing_columns}")
    
    df = df.copy()
    if "alt_loc" in df.columns:
        df.loc[df["alt_loc"] == "", "alt_loc"] = "A"
        df = df.loc[df["alt_loc"] == "A"]
    
    df = label_node_id(df, granularity)
    return df.loc[df["atom_name"] == granularity]

def plot_atom_point_cloud(processed_df: pd.DataFrame, output_path: str) -> None:
    """Creates a 3D scatter plot of atom coordinates. 

    Parameters:
    - processed_df (pd.DataFrame): Processed atomic data from remove_altloc. 
    - output_path(str): Path to save the HTML plot. 
    """
    atom_point_fig = px.scatter_3d(processed_df, x="x_coord", y="y_coord", z="z_coord", color="atom_name")
    atom_point_fig.update_traces(marker_size=4)
    atom_point_fig.update_layout(margin=dict(l=20, r=20, b=50, t=50))
    atom_point_fig.write_html(output_path)
    print(f"Point cloud HTML saved to: {output_path}")

def create_graphein_object(processed_df: pd.DataFrame, raw_pdb_df: pd.DataFrame, pdb_code: str, granularity: str = "CA") -> nx.Graph:
    """
    Initializes a Graphein object for 3D graphing, adds nodes, and returns it. This is required to use the Graphein package for visualization. 

    Parameters:
    - processed_df (pd.DataFrame): Processed DataFrame with alt loc removed from remove_altloc. 
    - raw_pdb_df (pd.DataFrame): Raw PDB data.
    - pdb_code (str): PDB code for the structure. Eg. icbn.pdb
    - granularity (str): Granularity for graph (default: "CA" for alpha carbon).

    Returns:
    - nx.Graph: Graphein graph object that will be required to plot backbone graphs. 
    """
    graphein_object = initialise_graph_with_metadata(
        protein_df=processed_df,
        raw_pdb_df=raw_pdb_df,
        pdb_code=pdb_code,
        granularity=granularity,
    )
    graphein_object = add_nodes_to_graph(graphein_object)
    return graphein_object

def plot_backbone_graph(graphein_object: nx.Graph, pdb_code: str = "Protein Structure", output_path: Optional[str] = None) -> None:
    """
    Connects backbone residues in a Graphein graph and visualizes the structure.

    Parameters:
    - graphein_object (nx.Graph): Protein graph.
    - pdb_code (str): PDB code or title for the graph.
    - output_path (str): Path to save the HTML plot.
    """
    for chain_id in graphein_object.graph["chain_ids"]:
        chain_residues = [
            (node, data) for node, data in graphein_object.nodes(data=True) if data["chain_id"] == chain_id
        ]
        for i in range(len(chain_residues) - 1):
            current_residue = chain_residues[i]
            next_residue = chain_residues[i + 1]

            same_chain = current_residue[1]["chain_id"] == next_residue[1]["chain_id"]
            adjacent_residues = abs(current_residue[1]["residue_number"] - next_residue[1]["residue_number"]) == 1

            if same_chain and adjacent_residues:
                if graphein_object.has_edge(current_residue[0], next_residue[0]):
                    graphein_object.edges[current_residue[0], next_residue[0]]["kind"].add("backbone_bond")
                else:
                    graphein_object.add_edge(current_residue[0], next_residue[0], kind={"backbone_bond"})

    plot = plotly_protein_structure_graph(
        graphein_object,
        colour_edges_by="kind",
        colour_nodes_by="seq_position",
        label_node_ids=False,
        plot_title=f"{pdb_code} Backbone Protein Graph",
        node_size_multiplier=1,
    )
    
    if output_path:
        plot.write_html(output_path)
        print(f"Backbone graph HTML saved to: {output_path}")

def open_html_in_browser(title: str, path: str) -> None:
    """
    Opens the HTML file from running the GUI in the system's default web browser.
    """
    try:
        if os.path.exists(path):
            print(f"Opening {title} at {path}")
            webbrowser.open(f"file://{os.path.abspath(path)}", new=2)
        else:
            error_msg = f"File not found: {path}"
            print(error_msg)
            messagebox.showerror("Error", error_msg)
    except Exception as e:
        error_msg = f"Error opening {title}: {str(e)}"
        print(error_msg)
        messagebox.showerror("Error", error_msg)

def process_pdb_file(pdb_path: str, frame: ttk.Frame) -> None:
    """Processes the PDB file and visualizes the atom point cloud and backbone graph."""
    if not os.path.exists(pdb_path):
        messagebox.showerror("Error", f"File not found: {pdb_path}")
        return

    try:
        output_dir = os.path.dirname(os.path.abspath(OUTPUT_FILES["point_cloud"]))
        os.makedirs(output_dir, exist_ok=True)
        
        parsed_df, header = parse_pdb(pdb_path)
        
        for file_type, file_path in OUTPUT_FILES.items():
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                except OSError as e:
                    print(f"Warning: Could not remove existing {file_type} file: {e}")
        
        point_cloud_path = os.path.abspath(OUTPUT_FILES["point_cloud"])
        backbone_graph_path = os.path.abspath(OUTPUT_FILES["backbone_graph"])
    
        plot_atom_point_cloud(parsed_df, point_cloud_path)
        altloc_removed_df = remove_altloc(parsed_df)
        
        pdb_code = os.path.basename(pdb_path)
        graphein_object = create_graphein_object(altloc_removed_df, parsed_df, pdb_code)
        plot_backbone_graph(graphein_object, pdb_code=pdb_code, output_path=backbone_graph_path)

        create_visualization_buttons(frame, point_cloud_path, backbone_graph_path)

    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")
        raise 

def select_pdb_file(frame: ttk.Frame) -> None:
    """
    Opens a file dialog to select a PDB file and processes it.
    Parameters:
    - frame (ttk.Frame): GUI frame to display the outputs.
    """
    pdb_path = filedialog.askopenfilename(filetypes=[("PDB Files", "*.pdb")])
    if pdb_path:
        process_pdb_file(pdb_path, frame)

def create_visualization_buttons(frame: ttk.Frame, point_cloud_path: str, backbone_graph_path: str) -> None:
    """
    Creates and packs visualization buttons.
    """
    tk.Button(
        frame, 
        text="View Point Cloud", 
        command=lambda: open_html_in_browser("Point Cloud", point_cloud_path)
    ).pack(pady=5)
    
    tk.Button(
        frame, 
        text="View 3D Backbone Graph", 
        command=lambda: open_html_in_browser("Backbone Graph", backbone_graph_path)
    ).pack(pady=5)

def create_main_window() -> Tuple[tk.Tk, ttk.Frame]:
    """
    Creates and configures the main application window.
    
    Returns:
        Tuple containing the root window and main frame
    """
    root = tk.Tk()
    root.title("Protein Structure Visualizer")
    root.geometry("400x300")
    
    main_frame = ttk.Frame(root, padding="10")
    main_frame.pack(fill=tk.BOTH, expand=True)
    
    select_button = ttk.Button(
        main_frame,
        text="Select PDB File",
        command=lambda: select_pdb_file(main_frame)
    )
    select_button.pack(pady=20)
    
    def on_closing():
        """Handles window closing event"""
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            print("Exiting...")
            try:
                root.quit()
                print("ProteinStructureVisualizer GUI stopping...")
            except Exception as e:
                print(f"Error during quit: {e}")
            finally:
                root.destroy()
                print("ProteinStructureVisualizer GUI has closed.")
    
    root.protocol("WM_DELETE_WINDOW", on_closing)
    
    return root, main_frame

def main():
    """Main entry point for the application."""
    try:
        root, main_frame = create_main_window()
        root.mainloop()
    except Exception as e:
        print(f"Error starting application: {e}")
        raise SystemExit(1)

if __name__ == "__main__":
    main()

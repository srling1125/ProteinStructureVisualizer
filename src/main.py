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
from typing import Any
from graphein.protein.visualisation import plotly_protein_structure_graph
from graphein.protein.graphs import label_node_id
from graphein.protein.graphs import initialise_graph_with_metadata
from graphein.protein.graphs import add_nodes_to_graph

def parse_pdb(pdb_path: Optional[str] = None, model_index: int = 1, parse_header: bool = True) -> pd.DataFrame:
    """
    Parse a PDB file into a DataFrame using PandasPdb. 

    Parameters:
    - pdb_path (str): Path to the PDB file.
    - model_index (int): Model index to extract (default: 1).
    - parse_header (bool): Whether to parse the PDB header (default: True). 

    Returns: 
    - pd.DataFrame: Concatenated DataFrame of ATOM and HETATM data. 
    - dict: Parsed PDB header or None if parse_header is False.
    """
    atomic_df = PandasPdb().read_pdb(pdb_path)
    if parse_header:
        header = parsePDBHeader(pdb_path)
    else:
        header = None
    atomic_df = atomic_df.get_model(model_index)
    if len(atomic_df.df["ATOM"]) == 0:
        raise ValueError(f"No model found for index: {model_index}")
    
    return pd.concat([atomic_df.df["ATOM"], atomic_df.df["HETATM"]]), header

def process_dataframe(df: pd.DataFrame, granularity="CA") -> pd.DataFrame:
    """
    Process the DataFrame for use with Graphein. 

    Parameters: 
    - df (pd.DataFrame): Input atomic data from parse_pdb. 
    - granularity (str): Granularity for graph (default is "CA" for alpha carbon).DS_Store

    Returns:
    - pd.DataFrame: Processed DataFrame.
    """
    if "alt_loc" in df.columns:
        df = df.copy()
        df["alt_loc"] = df["alt_loc"].replace("", "A")
        df = df.loc[df["alt_loc"] == "A"]
    df = label_node_id(df, granularity)
    df = df.loc[df["atom_name"] == granularity]
    return df 

def plot_atom_point_cloud(df: pd.DataFrame, output_path: str):
    """Create a 3D scatter plot of atom coordinates. 

    Parameters:
    - df (pd.DataFrame): Input atomic data.
    - output_path(str): Path to save the HTML plot.
    """
    fig = px.scatter_3d(df, x="x_coord", y="y_coord", z="z_coord", color="atom_name")
    fig.update_traces(marker_size=4)
    fig.update_layout(margin=dict(l=20, r=20, b=50, t=50))
    fig.write_html(output_path)
    print(f"Point cloud HTML saved to: {output_path}")

def create_graphein_object(processed_df: pd.DataFrame, raw_pdb_df: pd.DataFrame, pdb_code: str, granularity: str = "CA") -> Any:
    """
    Initializes a Graphein object for 3D graphing, adds nodes, and returns it.

    Parameters:
    - processed_df (pd.DataFrame): Processed DataFrame.
    - raw_pdb_df (pd.DataFrame): Raw PDB data.
    - pdb_code (str): PDB code for the structure.
    - granularity (str): Granularity for graph (default: "CA").

    Returns:
    - nx.Graph: Graphein graph object.
    """
    graphein_object = initialise_graph_with_metadata(
        protein_df=processed_df,
        raw_pdb_df=raw_pdb_df,
        pdb_code=pdb_code,
        granularity=granularity,
    )
    graphein_object = add_nodes_to_graph(graphein_object)
    return graphein_object

def plot_backbone_graph(graphein_object: nx.Graph, pdb_code="Protein Structure", output_path: str = None):
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

def open_html_in_browser(title, path):
    """
    Opens the HTML file in the system's default web browser.
    """
    if os.path.exists(path):
        print(f"Opening {title} at {path}")
        webbrowser.open(f"file://{path}", new=2)  
    else:
        print(f"File not found: {path}")
        messagebox.showerror("Error", f"File not found: {path}")

def process_pdb_file(pdb_path: str, frame: ttk.Frame):
    """Processes the PDB file and visualizes the atom point cloud and backbone graph.

    Parameters:
    - pdb_path (str): Path to the PDB file.
    - frame (ttk.Frame): GUI frame to display the outputs.
    - html_frame (HtmlFrame): Frame to embed the HTML visualizations.
    """
    try:
     # Step 1: Parse the PDB file
        atomic_df, header = parse_pdb(pdb_path)

        # Step 2: Create a 3D point cloud
        point_cloud_path = os.path.abspath("point_cloud.html")
        plot_atom_point_cloud(atomic_df, point_cloud_path)

        # Step 3: Process the PDB data
        processed_df = process_dataframe(atomic_df)

        # Step 4: Create Graphein graph
        graphein_object = create_graphein_object(processed_df, atomic_df, pdb_code=pdb_path.split("/")[-1])

        # Step 5: Plot Backbone Graph
        backbone_graph_path = os.path.abspath("backbone_graph.html")
        plot_backbone_graph(graphein_object, pdb_code=pdb_path.split("/")[-1], output_path=backbone_graph_path)

        tk.Button(frame, text="View Point Cloud", command=lambda: open_html_in_browser("Point Cloud", point_cloud_path)).pack(pady=5)
        tk.Button(frame, text="View 3D Backbone Graph", command=lambda: open_html_in_browser("Backbone Graph", backbone_graph_path)).pack(pady=5)

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def select_pdb_file(frame: ttk.Frame):
    """
    Opens a file dialog to select a PDB file and processes it.
    Parameters:
    - frame (ttk.Frame): GUI frame to display the outputs.
    """
    pdb_path = filedialog.askopenfilename(filetypes=[("PDB Files", "*.pdb")])
    if pdb_path:
        process_pdb_file(pdb_path, frame)
      
def on_closing():
    """
    Closes the GUI loop.
    """
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

root = tk.Tk()
root.title("Protein Structure Visualizer")
root.geometry("600x400")

root.protocol("WM_DELETE_WINDOW", on_closing)

frame = ttk.Frame(root, padding="10")
frame.pack(side=tk.LEFT, fill="y")

label = tk.Label(frame, text="Select a PDB file to visualize:")
label.pack(pady=5)

button = tk.Button(frame, text="Browse", command=lambda: select_pdb_file(frame))
button.pack(pady=5)

exit_button = tk.Button(frame, text="Exit", command=on_closing)
exit_button.pack(pady=5)

root.mainloop()

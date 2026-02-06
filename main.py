import customtkinter as ctk
from tkinter import filedialog
import pandas as pd
import pubchempy as pcp
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys
import os

# --- NEW DEPENDENCY IMPORTS ---
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcTPSA
except ImportError:
    # If RDKit is not installed, the app will exit gracefully with a message
    print("RDKit not found. Please run: pip install rdkit-pypi or conda install rdkit")
    sys.exit()

# --- CONFIGURATION ---
COMPOUNDS_PER_HEATMAP = 30  # Max number of compounds per plot for readability


# --- CORE BACKEND HELPER FUNCTIONS ---

def calculate_properties(smiles):
    """Calculates all required properties using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # RDKit property calculation
    props = {
        'MW': Descriptors.MolWt(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'LogP': Descriptors.MolLogP(mol),  # Use RDKit MolLogP
        'ROTB': Descriptors.NumRotatableBonds(mol),
        'PSA': CalcTPSA(mol),
        'NumRings': Descriptors.RingCount(mol),
        # Number of carbons and heteroatoms for Muegge
        'NumCarbons': mol.GetNumAtoms() - Descriptors.NumHeteroatoms(mol),
        'NumHeteroatoms': Descriptors.NumHeteroatoms(mol)
    }
    return props


def evaluate_compound_rules(props, rule_name):
    """
    Evaluates compound properties against the selected rule.
    Returns: (compliance_results, is_accepted)

    NOTE: Unicode symbols (<=, E) replaced with standard text (<=, in) for compatibility.
    """
    compliance_results = {}  # For the heatmap (keys are descriptive parameter names)
    is_accepted = False  # For overall filtering

    if rule_name == "Lipinski's Rule of Five":
        # Criteria: MW <= 500, HBA <= 10, HBD <= 5, LogP <= 5
        compliance_results = {
            'MW <= 500': props['MW'] <= 500,
            'HBA <= 10': props['HBA'] <= 10,
            'HBD <= 5': props['HBD'] <= 5,
            'LogP <= 5': props['LogP'] <= 5
        }
        # Accepted if compliance score >= 3 (1 violation max)
        is_accepted = sum(compliance_results.values()) >= 3

    elif rule_name == "Miles Congreve et al. RO3":
        # Criteria: ROTB <= 3, MW < 300, HBD <= 3, HBA <= 3, LogP <= 3, PSA <= 60
        compliance_results = {
            'MW < 300': props['MW'] < 300,
            'ROTB <= 3': props['ROTB'] <= 3,
            'HBA <= 3': props['HBA'] <= 3,
            'HBD <= 3': props['HBD'] <= 3,
            'LogP <= 3': props['LogP'] <= 3,
            'PSA <= 60': props['PSA'] <= 60
        }
        # Accepted if all 6 rules pass
        is_accepted = all(compliance_results.values())

    elif rule_name == "Muegge method":
        # Criteria: 9 rules using standard operators and 'in' for range
        compliance_results = {
            'MW in [200, 600]': (props['MW'] >= 200) and (props['MW'] <= 600),
            'LogP in [-2, 5]': (props['LogP'] >= -2) and (props['LogP'] <= 5),
            'PSA <= 150': props['PSA'] <= 150,
            'Rings <= 7': props['NumRings'] <= 7,
            'Carbons > 4': props['NumCarbons'] > 4,
            'Heteroatoms > 1': props['NumHeteroatoms'] > 1,
            'ROTB <= 15': props['ROTB'] <= 15,
            'HBA <= 10': props['HBA'] <= 10,
            'HBD <= 5': props['HBD'] <= 5
        }
        # Accepted if all 9 rules pass
        is_accepted = all(compliance_results.values())

    elif rule_name == "The rules of Veber et al.":
        # Criteria: ROTB <= 10, PSA <= 140
        compliance_results = {
            'ROTB <= 10': props['ROTB'] <= 10,
            'PSA <= 140': props['PSA'] <= 140
        }
        # Accepted if both rules pass
        is_accepted = all(compliance_results.values())

    return compliance_results, is_accepted


def filter_accepted_candidates(df_results_full):
    """Filters the full results DataFrame to accepted candidates based on the 'IsAccepted' column."""
    if df_results_full is not None and 'IsAccepted' in df_results_full.columns:
        return df_results_full[df_results_full['IsAccepted'] == True].copy()
    return pd.DataFrame()


# --- CUSTOMTKINTER APPLICATION CLASS ---

class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        # 1. Configuration and State
        self.title("DrugLikenZ")

        # Define window dimensions
        app_width = 1000
        app_height = 700

        # Get screen dimensions
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        # Calculate coordinates for centering
        x = (screen_width / 2) - (app_width / 2)
        y = (screen_height / 2) - (app_height / 2)

        # Set geometry (width x height + x_offset + y_offset)
        self.geometry(f"{app_width}x{app_height}+{int(x)}+{int(y)}")

        self.protocol("WM_DELETE_WINDOW", self.on_closing)

        try:
            # Use r-string for path consistency. Replace 'app_icon.ico' with your file name.
            # This only works for .ico files on Windows/macOS.
            self.iconbitmap(r'filter_list.ico')
        except Exception as e:
            pass

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.file_path = None

        # Define a larger font for navigation buttons/labels
        self.nav_font = ctk.CTkFont(size=14, weight="bold")

        # State variables
        self.df_results_full = None
        self.heatmap_chunks = []
        self.current_chunk_index = 0
        self.total_duplicates = 0
        self.canvas = None

        # 2. Sidebar Frame (Controls)
        self.sidebar_frame = ctk.CTkFrame(self, width=250, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")
        self.sidebar_frame.grid_columnconfigure(0, weight=1)

        self.logo_label = ctk.CTkLabel(self.sidebar_frame, text="DrugLikenZ",
                                       font=ctk.CTkFont(size=18, weight="bold"))
        self.logo_label.pack(pady=(20, 10))

        # Rule Selection
        self.rule_label = ctk.CTkLabel(self.sidebar_frame, text="Select Filtering Rule:",
                                       font=ctk.CTkFont(size=12, weight="bold"))
        self.rule_label.pack(pady=(10, 0))

        self.rule_options = [
            "Lipinski's Rule of Five",
            "Miles Congreve et al. RO3",
            "Muegge method",
            "The rules of Veber et al."
        ]
        self.rule_selection_combobox = ctk.CTkOptionMenu(self.sidebar_frame,
                                                         values=self.rule_options,
                                                         command=self.reset_ui_on_rule_change)
        self.rule_selection_combobox.set("Lipinski's Rule of Five")
        self.rule_selection_combobox.pack(padx=20, pady=5)

        # File Input
        self.file_button = ctk.CTkButton(self.sidebar_frame, text="Browse & Analyze File",
                                         command=self.open_file_dialog)
        self.file_button.pack(padx=20, pady=5)

        self.file_status = ctk.CTkLabel(self.sidebar_frame, text="No file selected.", wraplength=200)
        self.file_status.pack(pady=(0, 10))

        # Column Input
        self.column_label = ctk.CTkLabel(self.sidebar_frame, text="SMILES Column Name:")
        self.column_label.pack(pady=(10, 0))

        self.column_entry = ctk.CTkEntry(self.sidebar_frame, placeholder_text="SMILES", width=200)
        self.column_entry.insert(0, "SMILES")
        self.column_entry.pack(padx=20, pady=10)

        # Duplicates Counter
        self.duplicate_label = ctk.CTkLabel(self.sidebar_frame, text="Duplicates Removed: 0", text_color="orange")
        self.duplicate_label.pack(pady=(10, 5))

        # --- Plot Navigation Controls ---
        self.nav_frame = ctk.CTkFrame(self.sidebar_frame)
        self.nav_frame.pack(padx=20, pady=15)
        self.nav_frame.grid_columnconfigure((0, 2), weight=1)

        self.prev_button = ctk.CTkButton(self.nav_frame, text="◄", width=70,
                                         command=lambda: self.navigate_plot(-1), state="disabled",
                                         font=self.nav_font)  # Applied font
        self.prev_button.grid(row=0, column=0, padx=5)

        self.page_label = ctk.CTkLabel(self.nav_frame, text="Page 0/0", font=self.nav_font)  # Applied font
        self.page_label.grid(row=0, column=1, padx=5)

        self.next_button = ctk.CTkButton(self.nav_frame, text="►", width=70,
                                         command=lambda: self.navigate_plot(1), state="disabled",
                                         font=self.nav_font)  # Applied font
        self.next_button.grid(row=0, column=2, padx=5)

        # --- Export Controls ---
        self.export_frame = ctk.CTkFrame(self.sidebar_frame)
        self.export_frame.pack(padx=20, pady=15)

        self.export_data_button = ctk.CTkButton(self.export_frame, text="EXPORT ACCEPTED CSV", command=self.export_data,
                                                state="disabled")
        self.export_data_button.pack(padx=10, pady=5)

        self.export_plot_button = ctk.CTkButton(self.export_frame, text="Export Current Heatmap",
                                                command=lambda: self.export_heatmap("current"), state="disabled")
        self.export_plot_button.pack(padx=10, pady=5)

        self.export_all_plots_button = ctk.CTkButton(self.export_frame, text="Export ALL Heatmaps",
                                                     command=lambda: self.export_heatmap("all"), state="disabled")
        self.export_all_plots_button.pack(padx=10, pady=5)

        # 3. Main Content Frame (for Plotting)
        self.main_frame = ctk.CTkFrame(self)
        self.main_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        self.main_frame.grid_columnconfigure(0, weight=1)
        self.main_frame.grid_rowconfigure(0, weight=1)

        # --- MOVED STATUS MESSAGE HERE ---
        # Container for plot (top)
        self.plot_container = ctk.CTkFrame(self.main_frame, fg_color="transparent")
        self.plot_container.pack(fill="both", expand=True)

        self.main_label = ctk.CTkLabel(self.plot_container, text="Select a file and rule to begin analysis.",
                                       font=ctk.CTkFont(size=16))
        self.main_label.pack(pady=50)

        # Status Message (bottom of main frame)
        self.status_message = ctk.CTkLabel(self.main_frame, text="", wraplength=600, text_color="orange")
        self.status_message.pack(side="bottom", pady=10)

        self.cite_button = ctk.CTkButton(self.sidebar_frame,
                                         text="CITE THIS TOOL",
                                         fg_color="transparent",
                                         border_width=2,
                                         text_color=("gray10", "#DCE4EE"),
                                         command=self.show_citation)
        self.cite_button.pack(padx=20, pady=(20, 10), side="bottom")

    def center_child(self, child, width, height):
        """Calculates coordinates to center a child window relative to the App window."""
        self.update_idletasks()
        px = self.winfo_x()
        py = self.winfo_y()
        pw = self.winfo_width()
        ph = self.winfo_height()
        x = px + (pw // 2) - (width // 2)
        y = py + 30
        child.geometry(f"{width}x{height}+{int(x)}+{int(y)}")

    def on_closing(self):
        """
        Forcefully stops all background tasks to prevent 'application has been destroyed' errors.
        """
        # 1. Stop the mainloop immediately
        self.quit()

        # 2. Cancel any scheduled 'after' events (like DPI scaling checks)
        # This is the part that usually triggers the 'check_dpi_scaling' error
        for after_id in self.tk.eval('after info').split():
            self.after_cancel(after_id)

        # 3. Destroy the widgets
        self.destroy()

    # --- UI Action Methods ---

    def reset_ui_on_rule_change(self, selected_rule):
        """Called when a new rule is selected to force re-analysis if a file is already loaded."""
        if self.file_path and self.df_results_full is not None:
            self.run_screening()
        else:
            self.reset_ui_state()

    def open_file_dialog(self):
        """Opens a file dialog and runs screening immediately."""
        file_path = filedialog.askopenfilename(
            filetypes=[("Data files", "*.csv *.tsv *.txt")],
            title="Select SMILES File"
        )
        if file_path:
            self.file_path = file_path
            self.file_status.configure(text=f"Selected: {file_path.split('/')[-1]}", text_color="green")
            self.run_screening()
        else:
            self.file_status.configure(text="No file selected.", text_color="red")
            self.reset_ui_state()

    def run_screening(self):
        """Performs analysis: deduplication, screening, chunking, and initial plotting based on selected rule."""

        smiles_column_name = self.column_entry.get()
        selected_rule = self.rule_selection_combobox.get()
        self.status_message.configure(text=f"Processing with {selected_rule}...", text_color="blue")
        self.update_idletasks()

        try:
            # 1. Read and Deduplicate
            if self.file_path.endswith('.csv'):
                df_smiles = pd.read_csv(self.file_path)
            else:
                df_smiles = pd.read_csv(self.file_path, sep='\t')

            if smiles_column_name not in df_smiles.columns:
                raise ValueError(f"Column '{smiles_column_name}' not found.")

            initial_count = len(df_smiles)
            df_smiles.drop_duplicates(subset=[smiles_column_name], keep='first', inplace=True)
            self.total_duplicates = initial_count - len(df_smiles)
            self.duplicate_label.configure(text=f"Duplicates Removed: {self.total_duplicates}", text_color="teal")

            # 2. Screen Compounds
            all_results = []

            compliance_column_keys = None

            for index, row in df_smiles.iterrows():
                smiles = row[smiles_column_name]

                # a. Calculate ALL required properties
                props = calculate_properties(smiles)
                if props is None:
                    continue

                # b. Get Name from PubChem (optional)
                try:
                    # Get name from PubChem (maintains previous naming structure)
                    drug_can = pcp.get_compounds(smiles, 'smiles')[0]
                    compound_name = drug_can.synonyms[0] if drug_can.synonyms else smiles
                except IndexError:
                    compound_name = smiles

                # c. Evaluate Rules
                compliance_results, is_accepted = evaluate_compound_rules(props, selected_rule)

                # Store the descriptive column names for later use (they are consistent)
                if compliance_column_keys is None:
                    compliance_column_keys = list(compliance_results.keys())

                # Store the full results (all props + compliance)
                record = {
                    'Name': compound_name,
                    'SMILES': smiles,
                    'IsAccepted': is_accepted,
                    **props,
                    **compliance_results  # Compliance results for the heatmap
                }
                all_results.append(record)

            self.df_results_full = pd.DataFrame(all_results).set_index('Name')

            # 3. Filtering and Chunking
            accepted_count_df = self.df_results_full[self.df_results_full['IsAccepted'] == True]

            # Create the heatmap dataframe (only compliance columns)
            df_results_int = self.df_results_full[compliance_column_keys].astype(int)

            self.heatmap_chunks = [
                df_results_int.iloc[i:i + COMPOUNDS_PER_HEATMAP]
                for i in range(0, len(df_results_int), COMPOUNDS_PER_HEATMAP)
            ]
            self.current_chunk_index = 0

            if not self.heatmap_chunks:
                raise ValueError(f"No valid compounds found after screening using {selected_rule}.")

            # 4. Initial Plot and UI Enable
            if hasattr(self, 'main_label') and self.main_label.winfo_exists():
                self.main_label.destroy()

            self.show_plot(self.heatmap_chunks[0])
            self.update_navigation_controls()

            # Corrected status message
            self.status_message.configure(
                text=f"Screening Complete! Database reduced to {len(accepted_count_df)} drug-like molecules.",
                text_color="darkcyan")

            self.export_data_button.configure(state="normal")
            self.export_plot_button.configure(state="normal")
            self.export_all_plots_button.configure(state="normal")

        except Exception as e:
            self.status_message.configure(text=f"Error during screening: {e}", text_color="red")
            self.reset_ui_state()

    def show_plot(self, results):
        """Displays a single heatmap chunk, coloring failing compound names red, and rotating x-axis labels."""

        # Clear the container frame, but not the status message which is outside the container
        for widget in self.plot_container.winfo_children():
            widget.destroy()

        cmap = LinearSegmentedColormap.from_list('Custom', ['red', 'green'], N=2)

        num_compounds = results.shape[0]
        fig_height = max(5, num_compounds * 0.4)
        fig, ax = plt.subplots(figsize=(7, fig_height))

        # 1. Draw the Heatmap (cbar=False removes the scale)
        sns.heatmap(results, vmax=1, vmin=0, center=0.5, cmap=cmap,
                    linewidths=1, linecolor='lightgray', annot=False, ax=ax,
                    cbar=False)

        ax.set_title(f"{self.rule_selection_combobox.get()} Compliance", fontsize=10)
        ax.set_xlabel("", fontsize=8)  # Keep X-axis label for context
        ax.set_ylabel("")  # Remove Y-axis title completely

        # Rotate X-axis labels by 45 degrees
        ax.tick_params(axis='x', which='major', labelsize=8, rotation=45)

        # 2. Dynamic Y-axis Label Coloring Logic
        current_names = results.index.tolist()

        # Get the acceptance status for the current chunk from the main results
        accepted_status = self.df_results_full.loc[current_names, 'IsAccepted']

        for i, name in enumerate(current_names):
            is_accepted = accepted_status.iloc[i]

            # Failing compounds are where IsAccepted is False
            label_color = 'red' if not is_accepted else 'black'

            # Set the color for the tick label
            ax.get_yticklabels()[i].set_color(label_color)

        # 3. Final Matplotlib Configuration
        ax.tick_params(axis='y', which='major', labelsize=8, rotation=0)
        plt.tight_layout()

        # Embed the plot into the CONTAINER
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_container)
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(fill=ctk.BOTH, expand=True, padx=5, pady=5)
        self.canvas.draw()

    # --- (Navigation and Export methods remain unchanged) ---

    def navigate_plot(self, direction):
        """Handles navigation between heatmap chunks."""
        new_index = self.current_chunk_index + direction

        if 0 <= new_index < len(self.heatmap_chunks):
            self.current_chunk_index = new_index
            self.show_plot(self.heatmap_chunks[self.current_chunk_index])
            self.update_navigation_controls()

    def update_navigation_controls(self):
        """Updates page number and button states."""
        total_pages = len(self.heatmap_chunks)
        current_page = self.current_chunk_index + 1

        self.page_label.configure(text=f"Page {current_page}/{total_pages}")

        self.prev_button.configure(state="normal" if current_page > 1 else "disabled")
        self.next_button.configure(state="normal" if current_page < total_pages else "disabled")

    def export_data(self):
        """Exports accepted candidates to CSV, including SMILES and all properties."""
        if self.df_results_full is None:
            self.status_message.configure(text="Error: No results to export.", text_color="red")
            return

        # Use the filter function which now relies on the IsAccepted column
        df_accepted = filter_accepted_candidates(self.df_results_full)

        save_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            initialfile=f'{self.rule_selection_combobox.get().replace(" ", "_")}_accepted_candidates.csv',
            title="Save Accepted Candidates CSV"
        )
        if save_path:
            # We export all calculated properties and compliance status
            df_accepted.to_csv(save_path, index_label='Compound_Name')
            self.status_message.configure(text=f"Exported {len(df_accepted)} candidates to CSV.",
                                          text_color="green")
        else:
            self.status_message.configure(text="Export cancelled.", text_color="red")

    def export_heatmap(self, mode):
        """Exports the current heatmap or all heatmaps as PNGs."""
        if not self.heatmap_chunks:
            self.status_message.configure(text="Error: No heatmaps to export.", text_color="red")
            return

        current_rule_name = self.rule_selection_combobox.get().replace(" ", "_")

        if mode == "current":
            if hasattr(self, 'canvas') and self.canvas:
                save_path = filedialog.asksaveasfilename(
                    defaultextension=".png",
                    initialfile=f'{current_rule_name}_heatmap_page_{self.current_chunk_index + 1}.png',
                    title="Save Current Heatmap"
                )
                if save_path:
                    self.canvas.figure.savefig(save_path, bbox_inches='tight')
                    self.status_message.configure(text=f"Current heatmap saved to PNG.", text_color="green")
            else:
                self.status_message.configure(text="Error: Heatmap not currently displayed.", text_color="red")

        elif mode == "all":
            initial_dir = filedialog.askdirectory(title="Select Directory to Save All Heatmaps")
            if not initial_dir:
                self.status_message.configure(text="Export cancelled.", text_color="orange")
                return

            for i, chunk in enumerate(self.heatmap_chunks):
                current_names = chunk.index.tolist()
                accepted_status = self.df_results_full.loc[current_names, 'IsAccepted']

                # Re-generate the figure for saving
                fig, ax = plt.subplots(figsize=(7, max(5, chunk.shape[0] * 0.4)))
                sns.heatmap(chunk, vmax=1, vmin=0, center=0.5, cmap=sns.color_palette(['red', 'green'], 2),
                            linewidths=1, linecolor='lightgray', annot=False, ax=ax, cbar=False)

                ax.set_title(f"{self.rule_selection_combobox.get()} - Page {i + 1} of {len(self.heatmap_chunks)}",
                             fontsize=10)
                ax.set_ylabel("")  # Remove Y-axis title completely
                ax.tick_params(axis='x', which='major', labelsize=8, rotation=45)

                # Apply color to Y-axis labels
                for j, name in enumerate(current_names):
                    label_color = 'red' if not accepted_status.iloc[j] else 'black'
                    ax.get_yticklabels()[j].set_color(label_color)

                ax.tick_params(axis='y', which='major', labelsize=8, rotation=0)
                plt.tight_layout()

                output_file = f"{current_rule_name}_heatmap_page_{i + 1}_of_{len(self.heatmap_chunks)}.png"
                fig.savefig(os.path.join(initial_dir, output_file), bbox_inches='tight')
                plt.close(fig)

            self.status_message.configure(text=f"Exported {len(self.heatmap_chunks)} heatmaps.",
                                          text_color="darkgreen")

    def reset_ui_state(self):
        """Resets application state after an error or failed file load."""
        self.df_results_full = None
        self.heatmap_chunks = []
        self.current_chunk_index = 0
        self.total_duplicates = 0
        self.duplicate_label.configure(text="Duplicates Removed: 0", text_color="orange")
        self.page_label.configure(text="Page 0/0")
        self.prev_button.configure(state="disabled")
        self.next_button.configure(state="disabled")
        self.export_data_button.configure(state="disabled")
        self.export_plot_button.configure(state="disabled")
        self.export_all_plots_button.configure(state="disabled")

        # Restore the initial label in the container
        for widget in self.plot_container.winfo_children():
            widget.destroy()
        self.main_label = ctk.CTkLabel(self.plot_container, text="Select a file and rule to begin analysis.",
                                       font=ctk.CTkFont(size=14))
        self.main_label.pack(pady=50)

    def show_citation(self):
        """
        Opens a modal window displaying citations for TWO authors.
        """
        # --- 1. Data for Citation (Two Authors) ---
        author1_last = "Mahdaoui"
        author1_first = "Yazid"
        author1_initial = "Y."

        author2_last = "Zakkoumi"
        author2_first = "Hana"
        author2_initial = "H."

        year = "2026"
        title = "DrugLikenZ: An Open-Source Graphical User Interface for Batch Evaluation of Rule-of-Five Compliance"
        publisher = "Zenodo"
        doi = "https://doi.org/10.5281/zenodo.18500549"

        # --- 2. Format Logic ---
        def get_formatted_citations():
            return {
                "APA": f"{author1_last}, {author1_initial}, & {author2_last}, {author2_initial} ({year}). {title}. {publisher}. {doi}",
                "MLA": f"{author1_last}, {author1_first}, and {author2_first} {author2_last}. \"{title}.\" {publisher}, {year}. {doi}.",
                "Chicago": f"{author1_last}, {author1_first}, and {author2_first} {author2_last}. {title}. {publisher}, {year}. {doi}.",
                "Harvard": f"{author1_last}, {author1_initial} and {author2_last}, {author2_initial} {year}, '{title}', {publisher}, viewed [Date], <{doi}>.",
                "Vancouver": f"{author1_last} {author1_initial}, {author2_last} {author2_initial}. {title}. {publisher}; {year}. Available from: {doi}."
            }

        citations = get_formatted_citations()

        # --- 3. UI Setup ---
        citation_window = ctk.CTkToplevel(self)
        citation_window.title("Cite")

        # Centering logic for Citation Window
        self.center_child(citation_window, 650, 600)

        citation_window.attributes("-topmost", True)

        main_frame = ctk.CTkFrame(citation_window, fg_color="transparent")
        main_frame.pack(fill="both", expand=True, padx=20, pady=20)

        title_label = ctk.CTkLabel(main_frame, text="Cite", font=ctk.CTkFont(size=20, weight="bold"))
        title_label.pack(pady=(0, 20))

        for style, text in citations.items():
            style_frame = ctk.CTkFrame(main_frame, fg_color="transparent")
            style_frame.pack(fill="x", pady=(0, 15))
            style_label = ctk.CTkLabel(
                style_frame,
                text=style,
                width=80,
                anchor="w",
                text_color="gray",
                font=ctk.CTkFont(size=12, weight="bold")
            )
            style_label.pack(side="left", padx=(0, 10), anchor="n")

            citation_textbox = ctk.CTkTextbox(style_frame, height=50, wrap="word", activate_scrollbars=False)
            citation_textbox.insert("0.0", text)
            citation_textbox.configure(state="disabled")
            citation_textbox.pack(side="left", fill="x", expand=True)

        # --- 4. Export File Generators ---

        def generate_ris_content():
            """Generates RIS content for RefMan, RefWorks, Mendeley, Zotero."""
            return f"""TY  - COMP
AU  - {author1_last}, {author1_first}
AU  - {author2_last}, {author2_first}
TI  - {title}
PY  - {year}
PB  - {publisher}
UR  - {doi}
ER  - """

        def generate_enw_content():
            """Generates ENW content specifically for EndNote."""
            return f"""%0 Computer Program
%A {author1_last}, {author1_first}
%A {author2_last}, {author2_first}
%T {title}
%D {year}
%I {publisher}
%R {doi}
"""

        def export_citation(format_name):
            content = ""
            ext = ""

            # --- Logic for File Type Selection ---
            if format_name == "BibTeX":
                bibtex_entry = f"""@misc{{druglikenz_{year},
  author = {{{author1_last}, {author1_first} and {author2_last}, {author2_first}}},
  title = {{{title}}},
  year = {{{year}}},
  publisher = {{{publisher}}},
  doi = {{{doi}}}
}}"""
                self.clipboard_clear()
                self.clipboard_append(bibtex_entry)
                self.status_message.configure(text="BibTeX copied to clipboard!", text_color="green")
                return

            elif format_name == "EndNote":
                content = generate_enw_content()
                ext = ".enw"

            elif format_name in ["RefMan", "RefWorks"]:
                content = generate_ris_content()
                ext = ".ris"

            # --- Save File Dialog ---
            if content:
                citation_window.attributes("-topmost", False)

                save_path = filedialog.asksaveasfilename(
                    defaultextension=ext,
                    filetypes=[(f"{format_name} File", f"*{ext}"), ("All Files", "*.*")],
                    initialfile=f"DrugLikenZ_citation{ext}",
                    title=f"Export to {format_name}"
                )

                citation_window.attributes("-topmost", True)

                if save_path:
                    try:
                        with open(save_path, "w", encoding="utf-8") as f:
                            f.write(content)
                        self.status_message.configure(text=f"Saved citation to {save_path}", text_color="green")
                    except Exception as e:
                        self.status_message.configure(text=f"Error saving file: {e}", text_color="red")

        # --- 5. Export Buttons ---
        export_frame = ctk.CTkFrame(main_frame, fg_color="transparent")
        export_frame.pack(fill="x", pady=(20, 0), side="bottom")

        export_formats = ["BibTeX", "EndNote", "RefMan", "RefWorks"]
        for fmt in export_formats:
            btn = ctk.CTkButton(export_frame, text=fmt, fg_color="transparent",
                                text_color=("#3B8ED0", "#1F6AA5"), hover_color=("gray90", "gray25"),
                                anchor="w", width=60,
                                command=lambda f=fmt: export_citation(f))
            btn.pack(side="left", padx=(0, 15))


if __name__ == "__main__":
    # Ensure matplotlib doesn't conflict with Tkinter's threading
    plt.switch_backend('TkAgg')

    ctk.set_appearance_mode("System")
    ctk.set_default_color_theme("blue")

    app = App()

    try:
        app.mainloop()
    except Exception as e:
        print(f"Application encountered an error: {e}")
    finally:
        # Final safety check to ensure the process exits cleanly
        try:
            app.quit()
            app.destroy()
        except:
            pass
        # Force exit the script to stop background threads
        os._exit(0)
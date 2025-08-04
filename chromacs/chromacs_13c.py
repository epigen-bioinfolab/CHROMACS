import tkinter as tk
import tkinter.font as tkFont
from tkinter import filedialog, messagebox
from tkinter import ttk
from tkinter import simpledialog

import subprocess
import os
import threading
import glob
import re
import pandas as pd
import sys
from pkg_resources import resource_filename
import gzip
import shutil
import webbrowser

import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fpdf import FPDF
from PIL import Image
import time

import json
from collections import Counter
import networkx as nx

from collections import defaultdict

def load_custom_font(root):
    font_path = resource_filename("chromacs", "fonts/Roboto-Regular.ttf")
    if os.path.exists(font_path):
        try:
            roboto_font = tkFont.Font(family="Roboto", size=10)
            print(f"Successfully loaded custom font: {font_path}")
            return roboto_font
        except Exception as e:
            print(f"Error loading custom font: {e}")
            return None
    else:
        print(f"Custom font file not found at {font_path}. Using system default font.")
        return None


class ATACSeqPipeline:
    def __init__(self, root):
        self.root = root
        self.root.title("CHROMACS")
        self.root.geometry("950x800")

        self.roboto_font = load_custom_font(self.root)
        if self.roboto_font:
            style = ttk.Style()
            style.configure('.', font=self.roboto_font)
        else:
            style = ttk.Style()
            style.configure('.', font=("TkDefaultFont", 10))

        self._make_title_bar()

        self.main_frame = tk.Frame(self.root)
        self.main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))

        self.params = {
            "step1": {},
            "step2": {},
            "step3": {},
            "step4": {},
            "step5": {}
        }

        self.notebook = ttk.Notebook(self.main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.step1_frame = ttk.Frame(self.notebook)
        self.step2_frame = ttk.Frame(self.notebook)
        self.step3_frame = ttk.Frame(self.notebook)
        self.step4_frame = ttk.Frame(self.notebook)
        self.step5_frame = ttk.Frame(self.notebook)

        self.notebook.add(self.step1_frame, text="Step 1: Data Setup")
        self.notebook.add(self.step2_frame, text="Step 2: Quality Control")
        self.notebook.add(self.step3_frame, text="Step 3: Reference Selection")
        self.notebook.add(self.step4_frame, text="Step 4: Peak Caller Selection")
        self.notebook.add(self.step5_frame, text="Step 5: Analysis")

        for i in range(1, 5):
            self.notebook.tab(i, state="disabled")

        self.setup_step1_ui()
        self.setup_step2_ui()
        self.setup_step3_ui()
        self.setup_step4_ui()
        self.setup_step5_ui()

        self._setup_output_console()

        self.status_bar = tk.Label(self.main_frame, text="ChromAcS v0.1", bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X, padx=5, pady=(0, 5))

    def _make_title_bar(self):
        BAR_HEIGHT = 90
        self.title_bar = tk.Frame(self.root, bg="#800000", height=BAR_HEIGHT)
        self.title_bar.pack(fill=tk.X)
        self.title_bar.pack_propagate(False)


        left_frame = tk.Frame(self.title_bar, bg="#800000")
        left_frame.pack(side=tk.LEFT, padx=10, pady=0)
        left_frame.pack_propagate(True)

        chromacs_logo_path = resource_filename("chromacs", "assets/ChromAcS.png")
        if os.path.exists(chromacs_logo_path):
            try:
                img = tk.PhotoImage(file=chromacs_logo_path)
                max_h = BAR_HEIGHT - 20

                if img.height() > max_h:
                    factor = max(1, img.height() // max_h)
                    img = img.subsample(factor, factor)

                self.chromacs_logo_img = img
                logo_label = tk.Label(left_frame, image=img, bg="#800000", cursor="hand2")
                logo_label.image = img
                vertical_padding = (BAR_HEIGHT - img.height()) // 2
                logo_label.pack(side=tk.LEFT, pady=vertical_padding)

                def open_chromacs_site(event=None):
                    webbrowser.open("https://github.com/epigen-bioinfolab/CHROMACS/tree/main")

                logo_label.bind("<Button-1>", open_chromacs_site)

            except Exception as e:
                print(f"Could not load ChromAcS logo: {str(e)}")

        title_font = ("Roboto", 16, "bold") if self.roboto_font else ("Helvetica", 16, "bold")
        title_label = tk.Label(left_frame,
                               text="CHROMACS: Chromatin Accessibility Analysis Suite",
                               bg="#800000", fg="white", font=title_font)
        title_label.pack(side=tk.LEFT, padx=5, pady=(BAR_HEIGHT // 4, BAR_HEIGHT // 4))

        right_frame = tk.Frame(self.title_bar, bg="#800000")
        right_frame.pack(side=tk.RIGHT, padx=10, pady=0)
        right_frame.pack_propagate(True)

        lab_logo_path = resource_filename("chromacs", "assets/lab_logo.png")
        if os.path.exists(lab_logo_path):
            try:
                img = tk.PhotoImage(file=lab_logo_path)
                max_h = BAR_HEIGHT - 20

                if img.height() > max_h:
                    factor = max(1, img.height() // max_h)
                    img = img.subsample(factor, factor)

                self.lab_logo_img = img
                lab_logo_label = tk.Label(right_frame, image=img, bg="#800000", cursor="hand2")
                lab_logo_label.image = img
                vertical_padding = (BAR_HEIGHT - img.height()) // 2
                lab_logo_label.pack(side=tk.LEFT, pady=vertical_padding)

                def open_lab_website(event=None):
                    webbrowser.open("https://www.epigen-bioinfolab.com/")

                lab_logo_label.bind("<Button-1>", open_lab_website)

            except Exception as e:
                print(f"Could not load lab logo: {str(e)}")

        lab_link = tk.Label(right_frame,
                            text="Epigen-BioinfoLab",
                            font=("Roboto", 12, "bold") if self.roboto_font else ("Helvetica", 12, "bold"),
                            fg="white", bg="#800000", cursor="hand2")
        lab_link.pack(side=tk.LEFT, padx=5, pady=(BAR_HEIGHT // 4, BAR_HEIGHT // 4))
        lab_link.bind("<Button-1>", lambda e: webbrowser.open("https://www.epigen-bioinfolab.com/"))

    def _setup_output_console(self):
        console_frame = tk.Frame(self.main_frame, bd=1, relief=tk.SUNKEN)
        console_frame.pack(fill=tk.BOTH, expand=False, padx=5, pady=(5, 0))

        console_label = tk.Label(console_frame, text="Output Console", anchor=tk.W,
                                 font=("Roboto", 10, "bold") if self.roboto_font else ("Helvetica", 10, "bold"))
        console_label.pack(fill=tk.X, padx=5, pady=(2, 0))

        self.output_text = tk.Text(console_frame, wrap=tk.WORD, height=12,
                                   font=("Roboto", 9) if self.roboto_font else ("Courier", 9))
        self.output_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=2, pady=2)

        scrollbar = tk.Scrollbar(console_frame, command=self.output_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.output_text.config(yscrollcommand=scrollbar.set)
        self.output_text.config(state=tk.DISABLED)

        clear_btn = tk.Button(console_frame, text="Clear", command=self._clear_output,
                              width=6, font=("Roboto", 8) if self.roboto_font else ("Helvetica", 8))
        clear_btn.pack(side=tk.BOTTOM, padx=2, pady=2)

    def _clear_output(self):
        self.output_text.config(state=tk.NORMAL)
        self.output_text.delete(1.0, tk.END)
        self.output_text.config(state=tk.DISABLED)

    def get_timestamp(self):
        import datetime
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def clear_frame(self, frame):
        for widget in frame.winfo_children():
            widget.destroy()

# ======================= UI setup for the pipeline ====================================================================

    # step 1: selection of raw data ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def setup_step1_ui(self):
        # Output directory widgets
        tk.Label(self.step1_frame, text="Base Output Directory:", font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.output_dir = tk.Entry(self.step1_frame, width=50)
        self.output_dir.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.step1_frame, text="Browse", command=self.browse_output_dir, bg='grey').grid(row=0, column=2, padx=10, pady=5)
        tk.Button(self.step1_frame, text="Create", command=self.create_new_output_dir, bg='lightblue').grid(row=0, column=3, padx=10, pady=5)

        instruction_text_2 = (
            " [Assign the Base Output Directory where all the corresponding results will be saved. ]\n"
            " \"Create\" allows a new output directory; \"Browse\" to use an existing directory alreay created"
        )
        tk.Label(self.step1_frame, text=instruction_text_2, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=1, column=0, columnspan=4, padx=10, pady=5, sticky="w")

        tk.Label(self.step1_frame, text="Raw Data Directory:", font=(self.roboto_font, 10, 'bold')).grid(row=2, column=0, sticky="w", padx=10, pady=5)
        self.raw_data_dir = tk.Entry(self.step1_frame, width=50)
        self.raw_data_dir.grid(row=2, column=1, padx=10, pady=5)
        tk.Button(self.step1_frame, text="Browse", command=self.browse_raw_data, bg='grey').grid(row=2, column=2, padx=10, pady=5)
        instruction_text_1 = "[ If 'Select specific samples' is unchecked, ALL FASTQ files in the RAW DATA DIRECTORY will be processed. ]"
        tk.Label(self.step1_frame, text=instruction_text_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=3, column=0, columnspan=4, padx=10, pady=5, sticky="w")

        #selective sample selection
        self.select_samples_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            self.step1_frame,
            text="Select specific samples",
            variable=self.select_samples_var,
            command=self.toggle_sample_selection
        ).grid(row=4, column=0, columnspan=3, pady=5)

        self.sample_listbox = tk.Listbox(self.step1_frame, selectmode=tk.MULTIPLE, height=6, width=70)
        self.sample_listbox.grid(row=5, column=0, columnspan=3, padx=10, pady=5)
        self.sample_listbox.grid_remove()

        scrollbar = tk.Scrollbar(self.step1_frame, orient="vertical")
        scrollbar.config(command=self.sample_listbox.yview)
        scrollbar.grid(row=5, column=3, sticky="ns")
        self.sample_listbox.config(yscrollcommand=scrollbar.set)

        tk.Label(self.step1_frame, text="Number of Threads (Default- 8):", font=(self.roboto_font, 10, 'bold')).grid(
            row=6, column=0, sticky="w", padx=10, pady=5)
        self.threads_entry = tk.Entry(self.step1_frame, width=10)
        self.threads_entry.grid(row=6, column=1, padx=10, pady=5)
        self.threads_entry.insert(0, "8") # default = 8?

        tk.Button(self.step1_frame, text="Save & Next", command=self.save_step1_next, bg="yellow green").grid(
            row=7, column=0, columnspan=3, pady=10)

        # helper functions of step1_ui

    def browse_raw_data(self):
        directory = filedialog.askdirectory()
        if directory:
            self.raw_data_dir.delete(0, tk.END)
            self.raw_data_dir.insert(0, directory)
            if self.select_samples_var.get():
                self.populate_sample_list()

    def browse_output_dir(self):
        directory = filedialog.askdirectory()
        if directory:
            self.output_dir.delete(0, tk.END)
            self.output_dir.insert(0, directory)

    def create_new_output_dir(self):
        base_path = filedialog.askdirectory(title="Select Base Path to Create Directory")
        if base_path:
            new_dir_name = simpledialog.askstring("New Directory", "Enter new directory name:")
            if new_dir_name:
                new_dir_path = os.path.join(base_path, new_dir_name)
                try:
                    os.makedirs(new_dir_path, exist_ok=True)
                    self.output_dir.delete(0, tk.END)
                    self.output_dir.insert(0, new_dir_path)
                    self.update_output_gui(f"Created directory: {new_dir_path}\n")
                except Exception as e:
                    messagebox.showerror("Error", f"Could not create directory:\n{e}")

    def toggle_sample_selection(self):
        if self.select_samples_var.get():
            self.sample_listbox.grid()
            self.populate_sample_list()
        else:
            self.sample_listbox.grid_remove()

    def populate_sample_list(self):
        self.sample_listbox.delete(0, tk.END)
        directory = self.raw_data_dir.get()
        if directory:
            fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz")) + glob.glob(
                os.path.join(directory, "*.fq.gz"))
            for f in sorted(fastq_files):
                self.sample_listbox.insert(tk.END, os.path.basename(f))

    def save_step1_next(self):
        self.params["step1"]["raw_data_dir"] = self.raw_data_dir.get().strip()
        self.params["step1"]["base_output_dir"] = self.output_dir.get().strip()


        if self.select_samples_var.get():
            raw_data_dir = self.params["step1"]["raw_data_dir"]
            selected_files = [
                os.path.join(raw_data_dir, self.sample_listbox.get(i))
                for i in self.sample_listbox.curselection()
            ]
            self.params["step1"]["selected_samples"] = selected_files
        else:
            self.params["step1"]["selected_samples"] = None

        if not self.params["step1"]["raw_data_dir"] or not self.params["step1"]["base_output_dir"]:
            messagebox.showerror("Error", "Please fill in both directories.")
            return

        #base name to show in step2
        sample_names = set()
        raw_data_dir = self.params["step1"]["raw_data_dir"]
        selected_files = self.params["step1"].get("selected_samples", [])

        if not selected_files:
            selected_files = glob.glob(os.path.join(raw_data_dir, "*.fastq.gz")) + \
                             glob.glob(os.path.join(raw_data_dir, "*.fq.gz"))

        for filepath in selected_files:
            base = os.path.basename(filepath)

            match = re.match(r"^(.*?)(?:_R?[12](?:_001)?|_[12])\.f(?:ast)?q\.gz$", base)

            if match:
                sample_name = match.group(1)
            else:
                sample_name = re.sub(r'(_R?[12](?:_001)?|_[12])\.f(?:ast)?q\.gz$', '', base)

            sample_names.add(sample_name)

        self.params["step1"]["auto_sample_names"] = sorted(sample_names)

        #update thread as set by the user
        threads = self.threads_entry.get().strip()
        try:
            threads = int(threads)
            if threads < 1:
                raise ValueError("Thread count must be >= 1")
        except ValueError:
            threads = 8
            self.update_output_gui("Invalid thread count. Defaulting to 8.\n")

        self.params["threads"] = threads

        self.sample_listbox_step2.delete(0, tk.END)
        for name in self.params["step1"]["auto_sample_names"]:
            self.sample_listbox_step2.insert(tk.END, name)

        self.notebook.tab(1, state="normal")
        self.notebook.select(1)

    # step 2: trim galore and sample names +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def setup_step2_ui(self):
        tk.Label(self.step2_frame, text="Auto-Detected Sample Names:", font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, sticky="w", padx=10, pady=5)

        #read-only list
        self.sample_listbox_step2 = tk.Listbox(self.step2_frame, selectmode=tk.SINGLE, height=6, width=70)
        self.sample_listbox_step2.grid(row=0, column=1, padx=10, pady=5)

        instruction_1 = (
            "Please HOVER on this box and SCROLL to see more samples (if available)\n\n"
            "Samples auto-extracted from filenames.\n"
            "Format: SRRXXXXX (ignores _1/_2 suffixes).\n\n"
            "A list of allowed sample name types detected by ChromAcS-\n"
            "sample_1.fastq.gz / fq.gz;  sample_2.fastq.gz / fq.gz\n"
            "sample_R1.fastq.gz / fq.gz; sample_R2.fasq.gz / fq.gz\n"
            "sample_S1_L002_R1_001.fastq.gz/ fq.gz; sample_S1_L002_R2_001.fastq.gz/ fq.gz\n"
        )

        tk.Label(self.step2_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=1, column=1, columnspan=3, padx=10, pady=5, sticky="w")

        self.skip_trimming_var = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step2_frame, text="Skip Trimming (Use Raw Data for Alignment)",
                       variable=self.skip_trimming_var, font=(self.roboto_font, 10, 'bold')).grid(row=4, column=0, columnspan=2, padx=10, pady=5)

        instruction_3 = (
            "Checking this will skip the trimming done by Trim Galore.\n"
            "Raw data will be used directly for alignment.\n"
            "Please check the FastQC and MultiQC reports if you are unsure of this step. \n"
        )

        tk.Label(self.step2_frame, text=instruction_3, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=5, column=1, columnspan=3, padx=10, pady=5, sticky="w")

        tk.Button(self.step2_frame, text="Save & Next", command=self.save_step2_next, bg="yellow green").grid(row=7, column=1, pady=5)
        tk.Button(self.step2_frame, text="Back", command=lambda: self.notebook.select(0), bg='salmon').grid(row=7, column=0, pady=5)

        #helper functions of step2_ui

    def save_step2_next(self):

        self.params["step2"]["sample_names"] = self.params["step1"]["auto_sample_names"]
        self.params["step2"]["skip_trimming"] = self.skip_trimming_var.get()

        if not self.params["step2"]["sample_names"]:
            messagebox.showerror("Error", "No samples detected. Check Step 1.")
            return

        self.notebook.tab(2, state="normal")
        self.notebook.select(2)

    # Step 3: Alignment Settings +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def setup_step3_ui(self):

        self.genome_map = {
            # vertebrates
            "Human (GRCh38)": {
                "bowtie_ref": "GRCh38",
                "macs_size": "hs",
                "ensembl_species": "homo_sapiens",
                "ensembl_cap": "Homo_sapiens",
                "ensembl_assembly": "GRCh38"
            },
            "Mouse (GRCm39)": {
                "bowtie_ref": "GRCm39",
                "macs_size": "mm",
                "ensembl_species": "mus_musculus",
                "ensembl_cap": "Mus_musculus",
                "ensembl_assembly": "GRCm39"
            },
            "Rat (mRatBN7.2)": {
                "bowtie_ref": "mRatBN7.2",
                "macs_size": "2.9e9",
                "ensembl_species": "rattus_norvegicus",
                "ensembl_cap": "Rattus_norvegicus",
                "ensembl_assembly": "mRatBN7.2"
            },
            "Cow (ARS-UCD1.3)": {
                "bowtie_ref": "ARS-UCD1.3",
                "macs_size": "2.7e9",
                "ensembl_species": "bos_taurus",
                "ensembl_cap": "Bos_taurus",
                "ensembl_assembly": "ARS-UCD1.3"
            },
            "Pig (Sscrofa11.1)": {
                "bowtie_ref": "Sscrofa11.1",
                "macs_size": "2.5e9",
                "ensembl_species": "sus_scrofa",
                "ensembl_cap": "Sus_scrofa",
                "ensembl_assembly": "Sscrofa11.1"
            },
            "Chicken (GRCg7b)": {
                "bowtie_ref": "GRCg7b",
                "macs_size": "1.2e9",
                "ensembl_species": "gallus_gallus",
                "ensembl_cap": "Gallus_gallus",
                "ensembl_assembly": "bGalGal1.mat.broiler.GRCg7b"
            },
            "Chimpanzee (Pan_tro_3.0)": {
                "bowtie_ref": "Pan_tro_3.0",
                "macs_size": "3.3e9",
                "ensembl_species": "pan_troglodytes",
                "ensembl_cap": "Pan_troglodytes",
                "ensembl_assembly": "Pan_tro_3.0"
            },

            "Dog (ROS_Cfam_1.0)": {
                "bowtie_ref": "ROS_Cfam_1.0",
                "macs_size": "2.4e9",
                "ensembl_species": "canis_lupus_familiaris",
                "ensembl_cap": "Canis_lupus_familiaris",
                "ensembl_assembly": "ROS_Cfam_1.0"
            },
            "Goat (ARS1)": {
                "bowtie_ref": "ARS1",
                "macs_size": "2.9e9",
                "ensembl_species": "capra_hircus",
                "ensembl_cap": "Capra_hircus",
                "ensembl_assembly": "ARS1"
            },
            "Rabbit (OryCun2.0)": {
                "bowtie_ref": "OryCun2.0",
                "macs_size": "2.7e9",
                "ensembl_species": "oryctolagus_cuniculus",
                "ensembl_cap": "Oryctolagus_cuniculus",
                "ensembl_assembly": "OryCun2.0"
            },
            "Gorilla (gorGor4)": {
                "bowtie_ref": "gorGor4",
                "macs_size": "3.5e9",
                "ensembl_species": "gorilla_gorilla",
                "ensembl_cap": "Gorilla_gorilla",
                "ensembl_assembly": "gorGor4"
            },
            "Rhesus Macaque (Mmul_10)": {
                "bowtie_ref": "Mmul_10",
                "macs_size": "2.9e9",
                "ensembl_species": "macaca_mulatta",
                "ensembl_cap": "Macaca_mulatta",
                "ensembl_assembly": "Mmul_10"
            },

            # fish, amphibian
            "Zebrafish (GRCz11)": {
                "bowtie_ref": "GRCz11",
                "macs_size": "1.4e9",
                "ensembl_species": "danio_rerio",
                "ensembl_cap": "Danio_rerio",
                "ensembl_assembly": "GRCz11"
            },
            "Frog (UCB_Xtro_10.0)": {
                "bowtie_ref": "UCB_Xtro_10.0",
                "macs_size": "1.7e9",
                "ensembl_species": "xenopus_tropicalis",
                "ensembl_cap": "Xenopus_tropicalis",
                "ensembl_assembly": "UCB_Xtro_10.0"
            },
            "Indian Cobra (Nana_v5)": {
                "bowtie_ref": "Nana_v5",
                "macs_size": "1.8e9",
                "ensembl_species": "naja_naja",
                "ensembl_cap": "Naja_naja",
                "ensembl_assembly": "Nana_v5"
            },
            "Atlantic Salmon (Ssal_v3.1)": {
                "bowtie_ref": "Ssal_v3.1",
                "macs_size": "3.0e9",
                "ensembl_species": "salmo_salar",
                "ensembl_cap": "Salmo_salar",
                "ensembl_assembly": "Ssal_v3.1"
            },
            # invertebrates
            "C. elegans (WBcel235)": {
                "bowtie_ref": "WBcel235",
                "macs_size": "1.0e8",
                "ensembl_species": "caenorhabditis_elegans",
                "ensembl_cap": "Caenorhabditis_elegans",
                "ensembl_assembly": "WBcel235"
            },
            "Drosophila (BDGP6)": {
                "bowtie_ref": "BDGP6",
                "macs_size": "1.4e8",
                "ensembl_species": "drosophila_melanogaster",
                "ensembl_cap": "Drosophila_melanogaster",
                "ensembl_assembly": "BDGP6"
            },
        }

        tk.Label(self.step3_frame, text="Select Organism/Genome:", font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, sticky="w", padx=10,pady=5)
        self.genome_var = tk.StringVar()
        self.genome_menu = ttk.Combobox(self.step3_frame, textvariable=self.genome_var, values=list(self.genome_map.keys()), state="readonly")
        self.genome_menu.grid(row=0, column=1, padx=10, pady=5)

        self.ref_status = tk.Label(self.step3_frame, text="", fg="salmon")
        self.ref_status.grid(row=1, column=1, sticky="w", padx=10, pady=5)

        tk.Button(self.step3_frame, text="Check Reference", command=self.check_reference, bg="gray").grid(row=2, column=1, pady=5)
        tk.Button(self.step3_frame, text="Save & Next", command=self.save_step3_next, bg="yellow green").grid(row=3, column=1, pady=5)
        tk.Button(self.step3_frame, text="Back", command=lambda: self.notebook.select(1), bg="salmon").grid(row=3, column=0, pady=5)

    def check_reference(self):
        selected = self.genome_var.get()
        if not selected:
            messagebox.showerror("Error", "Please select a genome reference")
            return False

        genome_info = self.genome_map[selected]


        script_dir = os.path.dirname(os.path.abspath(__file__))
        ref_dir = os.path.join(script_dir, "ref_genome", genome_info["bowtie_ref"])

        bt2_base = os.path.join(ref_dir, genome_info["bowtie_ref"])
        fa_file = os.path.join(ref_dir, f"{genome_info['bowtie_ref']}.fa")
        fa_gz = f"{fa_file}.gz"


        index_files = [f"{bt2_base}.{ext}" for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]
        index_exists = all(os.path.exists(f) for f in index_files)

        if index_exists:
            self.ref_status.config(text="Reference index found!", fg="yellow green")
            return True
        else:
            status_messages = []
            if os.path.exists(fa_file):
                status_messages.append("FASTA file exists")
            elif os.path.exists(fa_gz):
                status_messages.append("Compressed FASTA exists")
            else:
                status_messages.append("No reference files found")

            status_text = "Reference status: " + ", ".join(status_messages) + ". Will be processed during pipeline run."
            self.ref_status.config(text=status_text, fg="orange")

            self.params["step3"].update({
                "genome_version": genome_info["bowtie_ref"],
                "ref_dir": ref_dir,
                "bt2_base": bt2_base,
                "ensembl_species": genome_info["ensembl_species"],
                "ensembl_cap": genome_info["ensembl_cap"],
                "ensembl_assembly": genome_info["ensembl_assembly"]
            })
            return False

    def save_step3_next(self):
        if not self.genome_var.get():
            messagebox.showerror("Error", "Select a genome reference.")
            return

        selected = self.genome_var.get()
        genome_info = self.genome_map[selected]

        script_dir = os.path.dirname(os.path.abspath(__file__))
        ref_dir = os.path.join(script_dir, "ref_genome", genome_info["bowtie_ref"])
        bt2_base = os.path.join(ref_dir, genome_info["bowtie_ref"])

        self.params["step3"].update({
            "genome_version": genome_info["bowtie_ref"],
            "macs3_genome_size": genome_info["macs_size"],
            "ensembl_species": genome_info["ensembl_species"],
            "ensembl_cap": genome_info["ensembl_cap"],
            "ensembl_assembly": genome_info["ensembl_assembly"],
            "ref_dir": ref_dir,
            "bt2_base": bt2_base
        })

        self.notebook.tab(3, state="normal")
        self.notebook.select(3)

    #step 4: peak caller selection +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def setup_step4_ui(self):
        tk.Label(self.step4_frame, text="Peak Calling", font=("Arial", 16, "bold")).grid(row=0, column=0,columnspan=2, pady=10)
        tk.Label(self.step4_frame, text="Select Peak Caller:", font=(self.roboto_font, 10, 'bold')).grid(row=1,column=0,sticky="w",padx=10,pady=5)
        self.peak_caller_choice = tk.StringVar(value="Genrich")
        tk.Radiobutton(self.step4_frame, text="Genrich", variable=self.peak_caller_choice, value="Genrich").grid(row=2,column=0,sticky="w",padx=10)
        tk.Radiobutton(self.step4_frame, text="MACS3", variable=self.peak_caller_choice, value="MACS3").grid(row=2,column=1,sticky="w",padx=10)
        tk.Button(self.step4_frame, text="Save & Next", command=self.save_step4_next, bg="yellow green").grid(row=3,column=2,columnspan=3,pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.notebook.select(2), bg='salmon').grid(row=3,column=0,pady=10)

    def save_step4_next(self):
        self.params["step4"]["peak_caller"] = self.peak_caller_choice.get().strip()
        #choice on peak caller
        if self.peak_caller_choice.get() == "Genrich":
            self.setup_step4a_genrich_ui()
        elif self.peak_caller_choice.get() == "MACS3":
            self.setup_step4b_macs3_ui()

    def setup_step4a_genrich_ui(self, samples=None):
        #genrich UI
        self.clear_frame(self.step4_frame)
        tk.Label(self.step4_frame, text="Genrich Peak Calling", font=("Arial", 16, "bold")).grid(
            row=0, column=0, columnspan=3, pady=10)
        tk.Label(self.step4_frame, text="Select Peak Type:", font=(self.roboto_font, 10, 'bold')).grid(row=1, column=0,sticky="w", padx=10, pady=5)
        self.peak_type_choice = tk.StringVar(value="merged")
        self.peak_type_choice.trace("w", lambda *args: self.toggle_merged_output_field())
        tk.Radiobutton(self.step4_frame, text="Merged", variable=self.peak_type_choice, value="merged").grid(row=1, column=1, sticky="w", padx=10)
        tk.Radiobutton(self.step4_frame, text="Unmerged", variable=self.peak_type_choice, value="unmerged").grid(row=1, column=2, sticky="w", padx=10)


        instruction_1 = (
            "If UNMERGED is checked, SEPARATE PEAK FILES are generated for each sample from the pool of multiple samples\n"
            "If MERGED is checked, ONE MERGED PEAK FILE is generated from the pool of multiple samples [suitable for replicates]"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=1, column=3, columnspan=4, padx=10, pady=5, sticky="w")


        step2_samples = self.params["step2"].get("sample_names", [])
        if not step2_samples:
            tk.Label(self.step4_frame, text="No sample names found. Please set them in Step 2.").grid(
                row=2, column=0, columnspan=3, padx=10, pady=5)
            return


        tk.Label(self.step4_frame, text="Test Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=3, column=0,sticky="w", padx=10, pady=5)
        self.treated_sample_var = tk.StringVar(value="Select a test sample")
        treated_dropdown = tk.OptionMenu(self.step4_frame, self.treated_sample_var, *step2_samples)
        treated_dropdown.grid(row=3, column=1, padx=10, pady=5)
        self.treated_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.treated_listbox.grid(row=3, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Test Sample", command=self.add_treated_sample, bg='gray').grid(
            row=3, column=3, padx=10, pady=5)


        instruction_1 = (
            "Select one by one: click on a sample and then click on \'Add Test Sample\'\n"
            "Repeat this for the next sample \n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=3, column=4, columnspan=4, padx=10, pady=5, sticky="w")


        tk.Label(self.step4_frame, text="Control Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=4, column=0,sticky="w", padx=10, pady=5)
        self.control_sample_var = tk.StringVar(value="Select a control sample")
        control_dropdown = tk.OptionMenu(self.step4_frame, self.control_sample_var, *step2_samples)
        control_dropdown.grid(row=4, column=1, padx=10, pady=5)
        self.control_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.control_listbox.grid(row=4, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Control Sample", command=self.add_control_sample, bg='gray').grid(row=4, column=3, padx=10, pady=5)

        instruction_1 = (
            "Note: You can skip this Control section, if there is none\n"
            "Please keep in mind that these are not biological controls, rather technical controls\n"
            "Do NOT put BASELINE condition as control here\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=4, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # blacklist files
        tk.Label(self.step4_frame, text="Select Blacklist Files:", font=(self.roboto_font, 10, 'bold')).grid(row=5,column=0,sticky="w",padx=10,pady=5)
        self.blacklist_listbox = tk.Listbox(self.step4_frame, selectmode=tk.MULTIPLE, height=5, width=50)
        self.blacklist_listbox.grid(row=5, column=1, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Files", command=self.browse_blacklist_files, bg='gray').grid(row=5,column=2,padx=10,pady=5)
        tk.Button(self.step4_frame, text="Remove Selected", command=self.remove_blacklist_files, bg='gray').grid(row=5,column=3,padx=10,pady=5)


        instruction_1 = (
            "Note: You can download blacklist from our Github page\n"
            "You can skip this part if you are unable to obtain a blacklist for your target organism\n"
            "Please ensure that you are using blacklists from the appropriate assembly versions\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=5, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        #exclude chromosomes
        tk.Label(self.step4_frame, text="Exclude Chromosomes:", font=(self.roboto_font, 10, 'bold')).grid(row=6,column=0,sticky="w",padx=10,pady=5)
        self.exclude_chr_m = tk.BooleanVar(value=False)
        self.exclude_chr_y = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step4_frame, text="Exclude MT (mitochondrial chromosomes)",variable=self.exclude_chr_m).grid(row=6, column=1,padx=10, pady=5)
        tk.Checkbutton(self.step4_frame, text="Exclude Y", variable=self.exclude_chr_y).grid(row=6, column=2,padx=10, pady=5)

        tk.Label(self.step4_frame, text="Custom Chromosomes for Exclusion").grid(row=6, column=3, sticky="w", padx=10,pady=5)
        self.custom_chr_entry = tk.Entry(self.step4_frame, width=30)
        self.custom_chr_entry.grid(row=6, column=4, padx=10, pady=5)

        # expand cutsites
        tk.Label(self.step4_frame, text='Expand cut-sites to _ bp:', font=(self.roboto_font, 10, 'bold')).grid(row=7,column=0,padx=10,pady=5,sticky='w')
        self.expand_cut_sites_var = tk.StringVar(value="100")
        self.expand_cut_combobox = ttk.Combobox(
            self.step4_frame,
            textvariable=self.expand_cut_sites_var,
            values=["100", "120", "150"]
        )
        self.expand_cut_combobox.grid(row=7, column=1, padx=10, pady=5)

        instruction_2 = (
            "For custom chromosome entry ensure to enter the correct chromosome naming format, like for Homo sapiens: 22\n"
            "For multiple entries they must be SINGLE WHITE-SPACED; like, 1 2 X 22\n"
            "Do NOT use prefix like 'chr'"
        )

        tk.Label(self.step4_frame, text=instruction_2, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=7, column=4, columnspan=6, padx=10, pady=5, sticky="w")

        tk.Label(self.step4_frame, text='Assign Max q-value (0 to 1):', font=(self.roboto_font, 10, 'bold')).grid(row=8,column=0,padx=10,pady=5)
        self.max_q_value_var = tk.StringVar(value="0.05")
        self.max_q_combobox = ttk.Combobox(
            self.step4_frame,
            textvariable=self.max_q_value_var,
            values=["0.05", "1"]
        )
        self.max_q_combobox.grid(row=8, column=1, padx=10, pady=5)

        self.merged_output_label = tk.Label(self.step4_frame, text="Merged Output Name:",font=(self.roboto_font, 10, 'bold'))
        self.merged_output_label.grid(row=9, column=0, sticky="w", padx=10, pady=5)
        self.merged_output_name_var = tk.StringVar()
        self.merged_output_entry = tk.Entry(self.step4_frame, textvariable=self.merged_output_name_var, width=30)
        self.merged_output_entry.grid(row=9, column=1, padx=10, pady=5)

        self.save_button = tk.Button(
            self.step4_frame,
            text="Save and Next",
            command=self.save_step4a_settings,
            bg="yellow green"
        )
        self.save_button.grid(row=12, column=0, columnspan=3, pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.navigate_back_to_step4_ui(), bg='salmon').grid(row=12,column=0,pady=10)

    def save_step4a_settings(self):
        self.params["step4"]["peak_type"] = self.peak_type_choice.get().strip()
        self.params["step4"]["treated_samples"] = [self.treated_listbox.get(i) for i in
                                                   self.treated_listbox.curselection()]
        self.params["step4"]["control_samples"] = [self.control_listbox.get(i) for i in
                                                   self.control_listbox.curselection()]
        self.params["step4"]["blacklist_files"] = [self.blacklist_listbox.get(i) for i in
                                                   range(self.blacklist_listbox.size())]


        exclude_chr_list = []
        if self.exclude_chr_m.get():
            exclude_chr_list.append("MT")
        if self.exclude_chr_y.get():
            exclude_chr_list.append("Y")

        custom_chr = self.custom_chr_entry.get().strip()
        if custom_chr:

            custom_list = custom_chr.split()
            exclude_chr_list.extend(custom_list)
        self.params["step4"]["exclude_chr"] = exclude_chr_list


        try:
            expand_cut = int(self.expand_cut_sites_var.get().strip())
            max_q = float(self.max_q_value_var.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Invalid numeric input for expand cut-sites or Q-value!")
            return

        self.params["step4"]["expand_cut_sites"] = expand_cut
        self.params["step4"]["max_q_value"] = max_q
        self.params["step4"]["merged_output_name"] = self.merged_output_name_var.get().strip()

        messagebox.showinfo("Step 4a Saved", "Genrich settings saved.")

        self.notebook.tab(4, state="normal")
        self.notebook.select(4)


    def toggle_merged_output_field(self):
        if self.peak_type_choice.get() == "merged":
            self.merged_output_label.grid(row=9, column=0, sticky="w", padx=10, pady=5)
            self.merged_output_entry.grid(row=9, column=1, padx=10, pady=5)

            self.save_button.grid(row=10, column=0, columnspan=3, pady=10)
        else:
            self.merged_output_label.grid_remove()
            self.merged_output_entry.grid_remove()
            self.save_button.grid(row=12, column=0, columnspan=3, pady=10)

    def add_treated_sample(self):
        sample = self.treated_sample_var.get()
        if sample and sample != "Select a test sample" and sample not in self.treated_listbox.get(0, tk.END):
            self.treated_listbox.insert(tk.END, sample)
            self.treated_listbox.select_set(tk.END)
        else:
            print("Test sample not added.")

    def add_control_sample(self):
        sample = self.control_sample_var.get()
        if sample and sample != "Select a control sample" and sample not in self.control_listbox.get(0, tk.END):
            self.control_listbox.insert(tk.END, sample)
            self.control_listbox.select_set(tk.END)
        else:
            print("Control sample not added.")

    def browse_blacklist_files(self):
        files = filedialog.askopenfilenames(
            title="Select Blacklist Files",
            filetypes=[("BED files", "*.bed"), ("All Files", "*.*")]
        )
        for file in files:
            if file not in self.blacklist_listbox.get(0, tk.END):
                self.blacklist_listbox.insert(tk.END, file)

    def remove_blacklist_files(self):
        selected_indices = self.blacklist_listbox.curselection()
        for idx in reversed(selected_indices):
            self.blacklist_listbox.delete(idx)

    def navigate_back_to_step4_ui(self):
        self.clear_frame(self.step4_frame)
        self.setup_step4_ui()

    def setup_step4b_macs3_ui(self, samples=None):
        #macs3 UI
        self.clear_frame(self.step4_frame)
        tk.Label(self.step4_frame, text="MACS3 Peak Calling", font=("Arial", 16, "bold")).grid(
            row=0, column=0, columnspan=3, pady=10)

        # uses sample names from step2
        if samples is None:
            step2_samples = self.params["step2"].get("sample_names", [])
        else:
            step2_samples = samples
        if not step2_samples:
            tk.Label(self.step4_frame, text="Enter Sample Names (separated by spaces):").grid(row=1, column=0, sticky="w", padx=10, pady=5)
            self.sample_input = tk.Entry(self.step4_frame, width=50)
            self.sample_input.grid(row=1, column=1, columnspan=2, padx=10, pady=5)
            tk.Button(self.step4_frame, text="Submit Samples", command=self.populate_samples).grid(row=1, column=3, padx=10, pady=5)
            return


        tk.Label(self.step4_frame, text="Test Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=2, column=0,sticky="w", padx=10,pady=5)
        self.treated_sample_var = tk.StringVar(value="Select a test sample")
        treated_dropdown = tk.OptionMenu(self.step4_frame, self.treated_sample_var, *step2_samples)

        treated_dropdown.grid(row=2, column=1, padx=10, pady=5)
        self.treated_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.treated_listbox.grid(row=2, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Test Sample", command=self.add_treated_sample, bg='gray').grid(row=2, column=3, padx=10, pady=5)

        instruction_1 = (
            "Select one by one: click on a sample and then click on \'Add Test Sample\'\n"
            "Repeat this for the next sample; you can ADD MULTIPLE SAMPLES at a time \n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=2, column=4, columnspan=4, padx=10, pady=5, sticky="w")


        tk.Label(self.step4_frame, text="Control Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=3, column=0,sticky="w",padx=10, pady=5)
        self.control_sample_var = tk.StringVar(value="Select a control sample")
        control_dropdown = tk.OptionMenu(self.step4_frame, self.control_sample_var, *step2_samples)
        control_dropdown.grid(row=3, column=1, padx=10, pady=5)

        self.control_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.control_listbox.grid(row=3, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Control Sample", command=self.add_control_sample, bg='gray').grid(row=3, column=3, padx=10, pady=5)

        instruction_1 = (
            "Note: You can skip this Control section, if there is none\n"
            "Please keep in mind that these are not biological controls, rather technical controls\n"
            "Do NOT put BASELINE condition as control here\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=3, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # blacklist file
        tk.Label(self.step4_frame, text="Select Blacklist Files:", font=(self.roboto_font, 10, 'bold')).grid(row=4,column=0,sticky="w",padx=10,pady=5)
        self.blacklist_listbox = tk.Listbox(self.step4_frame, selectmode=tk.MULTIPLE, height=5, width=50)
        self.blacklist_listbox.grid(row=4, column=1, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Files", command=self.browse_blacklist_files, bg='gray').grid(row=4,column=2,padx=10,pady=5)
        tk.Button(self.step4_frame, text="Remove Selected", command=self.remove_blacklist_files, bg='gray').grid(row=4,column=3,padx=10,pady=5)


        instruction_1 = (
            "Note: You can download blacklist from our Github page\n"
            "You can skip this part if you are unable to obtain a blacklist for your target organism\n"
            "Please ensure that you are using blacklists from the appropriate assembly versions\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=4, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        tk.Label(self.step4_frame, text="Exclude Chromosomes:", font=(self.roboto_font, 10, 'bold')).grid(row=5,column=0,sticky="w",padx=10,pady=5)
        self.exclude_chr_m = tk.BooleanVar(value=False)
        self.exclude_chr_y = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step4_frame, text="Exclude MT (Mitochondrial chromosomes)", variable=self.exclude_chr_m).grid(row=5, column=1,padx=10, pady=5)
        tk.Checkbutton(self.step4_frame, text="Exclude Y", variable=self.exclude_chr_y).grid(row=5, column=2,padx=10, pady=5)
        tk.Label(self.step4_frame, text="Custom Chromosomes for Exclusion").grid(row=5, column=3, sticky="w", padx=10, pady=5)
        self.custom_chr_entry = tk.Entry(self.step4_frame, width=30)
        self.custom_chr_entry.grid(row=5, column=4, padx=10, pady=5)


        tk.Label(self.step4_frame, text='Assign Max q-value (0 to 1):', font=(self.roboto_font, 10, 'bold')).grid(row=6, column=0,padx=10,pady=5,sticky="w")
        self.max_q_value_var = tk.StringVar(value="0.05")
        self.max_q_combobox = ttk.Combobox(
            self.step4_frame,
            textvariable=self.max_q_value_var,
            values=["0.05", "0.01"]
        )
        self.max_q_combobox.grid(row=6, column=1, padx=10, pady=5)

        instruction_2 = (
            "For custom chromosome entry ensure to enter the correct chromosome naming format, like for Homo sapiens: 22\n"
            "For multiple entries they must be SINGLE WHITE-SPACED; like, 1 2 X 22\n"
            "Do NOT use prefix like 'chr'"
        )

        tk.Label(self.step4_frame, text=instruction_2, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=6, column=4, columnspan=6, padx=10, pady=5, sticky="w")

        tk.Label(self.step4_frame, text="Genome Size (auto-filled):").grid(row=7, column=0, sticky="w", padx=10, pady=5)
        self.macs3_genome_entry = tk.Entry(self.step4_frame, width=30, state='readonly')
        self.macs3_genome_entry.grid(row=7, column=1, padx=10, pady=5)

        # auto filled from step3
        genome_size = self.params["step3"].get("macs3_genome_size", "hs")
        self.macs3_genome_entry.config(state='normal')
        self.macs3_genome_entry.insert(0, genome_size)
        self.macs3_genome_entry.config(state='readonly')
        tk.Button(self.step4_frame, text="Save and Next", command=self.save_step4b_settings, bg="yellow green").grid(row=9, column=0, columnspan=3, pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.navigate_back_to_step4_ui(), bg='salmon').grid(row=9,column=0,pady=10)

    def populate_samples(self):
        sample_input = self.sample_input.get()
        samples = sample_input.split()
        self.setup_step4b_macs3_ui(samples)

    def save_step4b_settings(self):
        self.params["step4"]["test_samples"] = [self.treated_listbox.get(i) for i in
                                                   self.treated_listbox.curselection()]
        self.params["step4"]["control_samples"] = [self.control_listbox.get(i) for i in
                                                   self.control_listbox.curselection()]
        self.params["step4"]["blacklist_files"] = [self.blacklist_listbox.get(i) for i in
                                                   range(self.blacklist_listbox.size())]
        exclude_chr_list = []
        if self.exclude_chr_m.get():
            exclude_chr_list.append("MT")
        if self.exclude_chr_y.get():
            exclude_chr_list.append("Y")
        custom_chr = self.custom_chr_entry.get().strip()
        if custom_chr:
            custom_list = custom_chr.split()
            exclude_chr_list.extend(custom_list)
        self.params["step4"]["exclude_chr"] = exclude_chr_list

        try:
            max_q = float(self.max_q_value_var.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Invalid Q-value!")
            return

        self.params["step4"]["max_q_value"] = max_q

        self.params["step4"]["genome"] = self.params["step3"]["macs3_genome_size"]

        messagebox.showinfo("Step 4b Saved", "MACS3 settings saved.")

        self.notebook.tab(4, state="normal")
        self.notebook.select(4)

    def setup_step5_ui(self):
        frame = self.step5_frame

        for i in range(3):
            frame.grid_columnconfigure(i, weight=1)
        frame.grid_rowconfigure(0, weight=1)

        instruction_1 = (
            "Click 'Run Pipeline' to generate peak files with annotation.\n"
            "If you are running ChromAcS for the first time, it may take longer due to reference genome setup and indexing.\n\n"
            "- FastQC, MultiQC, Trimmed data, Aligned BAM files, and Coverage files will be stored under the base output directory.\n"
            "- Peak files will be saved in the 'peak_files' folder.\n"
            "- Annotated peaks will be in 'Annotated_Peaks' inside 'peak_files'."
        )
        tk.Label(
            frame, text=instruction_1, wraplength=800, justify=tk.LEFT,
            anchor="center", font=(self.roboto_font, 10, 'bold')
        ).grid(row=0, column=1, padx=20, pady=10, sticky="ew")

        button_frame = tk.Frame(frame)
        button_frame.grid(row=1, column=1, pady=10)

        tk.Button(button_frame, text="Back", command=lambda: self.notebook.select(3), bg='salmon'
                  ).pack(side=tk.LEFT, padx=20)
        tk.Button(button_frame, text="Run Pipeline", command=self.run_pipeline, bg="yellow green"
                  ).pack(side=tk.LEFT, padx=20)

        #diffbind module
        diffbind_info = (
            "After peak calling is complete, you can run DiffBind to perform differential binding analysis "
            "based on experimental conditions."
        )
        tk.Label(
            frame, text=diffbind_info, wraplength=800, justify=tk.LEFT,
            anchor="center", font=(self.roboto_font, 9, 'italic')
        ).grid(row=2, column=1, padx=20, pady=(20, 5), sticky="ew")

        diffbind_btn_frame = tk.Frame(frame)
        diffbind_btn_frame.grid(row=3, column=1, pady=10)
        self.diffbind_btn = tk.Button(
            diffbind_btn_frame, text="Run DiffBind Analysis", bg='orange',
            command=self.launch_diffbind_config, state=tk.NORMAL
        )
        self.diffbind_btn.pack()

        #noiseq module
        noisq_info = (
            "If you do not have biological replicates, you can use NOISeq instead of DiffBind.\n"
            "NOISeq performs a simulation-based differential analysis on a peak count matrix."
        )
        tk.Label(
            frame, text=noisq_info, wraplength=800, justify=tk.LEFT,
            anchor="center", font=(self.roboto_font, 9, 'italic')
        ).grid(row=4, column=1, padx=20, pady=(10, 5), sticky="ew")

        # NOISeq button (centered below)
        noisq_btn_frame = tk.Frame(frame)
        noisq_btn_frame.grid(row=5, column=1, pady=(5, 20))
        self.noisq_btn = tk.Button(
            noisq_btn_frame, text="Run NOISeq Analysis (No Replicates)", bg='light blue',
            command=self.launch_noisq_config
        )
        self.noisq_btn.pack()

        for i in range(5):
            frame.grid_rowconfigure(i, weight=1)

        # motif enrichment module
        motif_info = (
            "After differential peak analysis, you can run Motif Enrichment to identify enriched TF motifs "
            "You need to download suitable .meme file from JASPAR manually"
        )
        tk.Label(
            frame, text=motif_info, wraplength=800, justify=tk.LEFT,
            anchor="center", font=(self.roboto_font, 9, 'italic')
        ).grid(row=6, column=1, padx=20, pady=(10, 5), sticky="ew")

        motif_btn_frame = tk.Frame(frame)
        motif_btn_frame.grid(row=7, column=1, pady=(5, 20))
        self.motif_btn = tk.Button(
            motif_btn_frame, text="Run Motif Enrichment", bg='medium orchid',
            command=self.launch_motif_module
        )
        self.motif_btn.pack()

        for i in range(8):
            frame.grid_rowconfigure(i, weight=1)


# ======================= Final Run Pipeline ===========================================================================

    def run_pipeline(self):
        base_out = self.params["step1"].get("base_output_dir")
        if not base_out:
            messagebox.showerror("Error", "Base output directory not specified in Step 1!")
            return

        fastqc_raw = os.path.join(base_out, "fastqc_raw")
        multiqc_raw = os.path.join(base_out, "multiqc_raw")
        trimmed_data = os.path.join(base_out, "trimmed_data")
        bam_output = os.path.join(base_out, "bam_output")
        fastqc_trimmed = os.path.join(base_out, "fastqc_trimmed")
        multiqc_trimmed = os.path.join(base_out, "multiqc_trimmed")
        normalized_coverage = os.path.join(base_out, "normalized_coverage")
        peak_files = os.path.join(base_out, "peak_files")
        os.makedirs(fastqc_raw, exist_ok=True)
        os.makedirs(multiqc_raw, exist_ok=True)
        os.makedirs(trimmed_data, exist_ok=True)
        os.makedirs(bam_output, exist_ok=True)
        os.makedirs(fastqc_trimmed, exist_ok=True)
        os.makedirs(multiqc_trimmed, exist_ok=True)
        os.makedirs(normalized_coverage, exist_ok=True)
        os.makedirs(peak_files, exist_ok=True)

        if not self.params["step5"].get("peak_files_dir"):
            self.params["step5"]["peak_files_dir"] = peak_files

        threading.Thread(target=self.execute_pipeline_sequentially, daemon=True).start()

    def execute_pipeline_sequentially(self):
        base_out = self.params["step1"]["base_output_dir"]
        raw_data = self.params["step1"]["raw_data_dir"]
        samples = self.params["step2"].get("sample_names", [])

        fastqc_raw = os.path.join(base_out, "fastqc_raw")
        multiqc_raw = os.path.join(base_out, "multiqc_raw")
        trimmed_data = os.path.join(base_out, "trimmed_data")
        bam_output = os.path.join(base_out, "bam_output")
        fastqc_trimmed = os.path.join(base_out, "fastqc_trimmed")
        multiqc_trimmed = os.path.join(base_out, "multiqc_trimmed")
        normalized_coverage = os.path.join(base_out, "normalized_coverage")
        peak_files = os.path.join(base_out, "peak_files")
        threads = self.params.get("threads", 8)

        try:

            selected_samples = self.params["step1"].get("selected_samples")

            if selected_samples:
                selected_names = {
                    os.path.basename(f).split("_")[0] for f in selected_samples
                }
                samples = [s for s in self.params["step1"]["auto_sample_names"] if s in selected_names]
            else:
                samples = self.params["step1"]["auto_sample_names"]

            # step1: fastqc on raw data
            self.update_output_gui("Checking Step 1: FastQC raw data...\n")
            raw_files = []
            if self.params["step1"].get("selected_samples"):
                raw_files = self.params["step1"]["selected_samples"]
                samples = self.params["step1"]["auto_sample_names"]
            else:
                raw_files = glob.glob(os.path.join(raw_data, "*.fastq.gz")) + glob.glob(
                    os.path.join(raw_data, "*.fq.gz"))

            all_raw_reports_exist = all(
                os.path.exists(os.path.join(fastqc_raw, os.path.basename(f).replace('.fastq.gz', '_fastqc.html').replace('.fq.gz', '_fastqc.html')))
                for f in raw_files
            )

            if not all_raw_reports_exist:
                self.update_output_gui("Running Step 1: FastQC on raw data...\n")
                files_to_process = " ".join(raw_files)
                cmd = f"fastqc -t {threads} {files_to_process} -o {fastqc_raw}"
                if not self.run_blocking_command(cmd):
                    return
            else:
                self.update_output_gui("Skipping Step 1: FastQC raw data (reports exist)\n")

            # step2: multiqc on raw data
            self.update_output_gui("Checking Step 2: MultiQC raw data...\n")
            multiqc_report = os.path.join(multiqc_raw, "multiqc_report.html")
            if not os.path.exists(multiqc_report):
                self.update_output_gui("Running Step 2: MultiQC on raw data...\n")
                cmd = f"multiqc {fastqc_raw} -o {multiqc_raw}"
                if not self.run_blocking_command(cmd):
                    return
            else:
                self.update_output_gui("Skipping Step 2: MultiQC raw data (report exists)\n")

            selected_samples = self.params["step1"].get("selected_samples")
            if not selected_samples:
                raw_data_dir = self.params["step1"]["raw_data_dir"]
                selected_samples = glob.glob(os.path.join(raw_data_dir, "*.fastq.gz")) + \
                                   glob.glob(os.path.join(raw_data_dir, "*.fq.gz"))

            # step3: trimming
            if not self.params["step2"].get("skip_trimming", False):
                self.update_output_gui("Checking Step 3: Trim Galore...\n")
                for sample in samples:
                    output_r1 = os.path.join(trimmed_data, f"{sample}_val_1.fq.gz")
                    output_r2 = os.path.join(trimmed_data, f"{sample}_val_2.fq.gz")

                    if not (os.path.exists(output_r1) and os.path.exists(output_r2)):
                        sample_files = [f for f in selected_samples if sample in os.path.basename(f)]

                        r1_list = [f for f in sample_files if
                                   re.search(r"_R?1(?:_001)?\.f(ast)?q\.gz$", os.path.basename(f))]
                        r2_list = [f for f in sample_files if
                                   re.search(r"_R?2(?:_001)?\.f(ast)?q\.gz$", os.path.basename(f))]

                        if not r1_list or not r2_list:
                            self.update_output_gui(f"Warning: Could not find both pairs for sample {sample}\n")
                            continue

                        r1 = r1_list[0]
                        r2 = r2_list[0]

                        self.update_output_gui(f"Running Trim Galore for {sample}...\n")
                        cmd = (f"trim_galore --basename {sample} --gzip  --cores {threads} --paired "
                               f"{r1} {r2} --output_dir {trimmed_data}")

                        if not self.run_blocking_command(cmd, show_output=False):
                            return
                    else:
                        self.update_output_gui(f"Skipping Trim Galore for {sample} (trimmed files exist)\n")
            else:
                self.update_output_gui("Skipping Step 3: Trim Galore (User chose to use raw data)\n")

            # step4: fastqc on trimmed data
            if not self.params["step2"].get("skip_trimming", False):
                self.update_output_gui("Checking Step 4: FastQC trimmed data...\n")
                trimmed_files = []
                for sample in samples:
                    trimmed_files.append(os.path.join(trimmed_data, f"{sample}_val_1.fq.gz"))
                    trimmed_files.append(os.path.join(trimmed_data, f"{sample}_val_2.fq.gz"))

                all_trimmed_reports_exist = all(
                    os.path.exists(os.path.join(fastqc_trimmed, os.path.basename(f).replace('.fq.gz', '_fastqc.html')))
                    for f in trimmed_files
                )
                if not all_trimmed_reports_exist:
                    self.update_output_gui("Running Step 4: FastQC on trimmed data...\n")
                    cmd = f"fastqc -t {threads} {trimmed_data}/*.fq.gz -o {fastqc_trimmed}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui("Skipping Step 4: FastQC trimmed data (reports exist)\n")
            else:
                self.update_output_gui(
                    "Skipping Steps 4 & 5: FastQC and MultiQC on trimmed data (trimming was skipped).\n")

            # step5: multiqc on trimmed data
            if not self.params["step2"].get("skip_trimming", False):
                self.update_output_gui("Checking Step 5: MultiQC trimmed data...\n")
                multiqc_report = os.path.join(multiqc_trimmed, "multiqc_report.html")
                if not os.path.exists(multiqc_report):
                    self.update_output_gui("Running Step 5: MultiQC on trimmed data...\n")
                    cmd = f"multiqc {fastqc_trimmed} -o {multiqc_trimmed}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui("Skipping Step 5: MultiQC trimmed data (report exists)\n")

            # ====== ensembl reference genome and building etc. ===============
            ensembl_release = "111"
            sp = self.params["step3"]["ensembl_species"]
            cap = self.params["step3"]["ensembl_cap"]
            assembly = self.params["step3"]["ensembl_assembly"]
            ref_dir = self.params["step3"]["ref_dir"]
            bt2_base = self.params["step3"]["bt2_base"]
            genome_version = self.params["step3"]["genome_version"]
            fa_file = os.path.join(ref_dir, f"{genome_version}.fa")
            fa_gz = f"{fa_file}.gz"
            index_extensions = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

            if not all(os.path.exists(f"{bt2_base}.{ext}") for ext in index_extensions):
                self.update_output_gui(f"Processing reference genome: {genome_version}\n")

                os.makedirs(ref_dir, exist_ok=True)


                if not os.path.exists(fa_file):
                    if os.path.exists(fa_gz):
                        self.update_output_gui("Found compressed FASTA, decompressing...\n")
                        subprocess.run(f"gzip -d {fa_gz}", shell=True, check=True)
                    else:
                        self.update_output_gui("Downloading reference genome with aria2c...\n")
                        url = (f"ftp://ftp.ensembl.org/pub/release-{ensembl_release}/fasta/"
                               f"{sp}/dna/{cap}.{assembly}.dna.toplevel.fa.gz")

                        download_cmd = (
                            f"aria2c -c -x 16 -s 16 -d {ref_dir} "
                            f"-o {genome_version}.fa.gz {url}"
                        )
                        try:
                            subprocess.run(download_cmd, shell=True, check=True)
                            subprocess.run(f"gzip -d {fa_gz}", shell=True, check=True)
                        except subprocess.CalledProcessError as e:
                            self.show_error_gui(f"Download failed: {str(e)}")
                            return


                if os.path.exists(fa_file):
                    self.update_output_gui("Building Bowtie2 index...\n")
                    try:
                        subprocess.run(
                            f"bowtie2-build --threads {threads} {fa_file} {bt2_base}",
                            shell=True,
                            check=True
                        )
                    except subprocess.CalledProcessError as e:
                        self.show_error_gui(f"Index build failed: {str(e)}")
                        return
                else:
                    self.show_error_gui("FASTA file missing after download/decompression")
                    return
            else:
                self.update_output_gui("Using existing reference index\n")

            # ======== download GTF =========================================
            gtf_url = (
                f"ftp://ftp.ensembl.org/pub/release-{ensembl_release}/gtf/"
                f"{sp}/{cap}.{assembly}.{ensembl_release}.gtf.gz"
            )
            gtf_gz = os.path.join(ref_dir, f"{genome_version}.gtf.gz")
            gtf_file = os.path.join(ref_dir, f"{genome_version}.gtf")

            if not os.path.exists(gtf_file):
                if os.path.exists(gtf_gz):
                    self.update_output_gui("Found compressed GTF, decompressing...\n")
                    subprocess.run(f"gzip -d {gtf_gz}", shell=True, check=True)
                else:
                    self.update_output_gui("Downloading GTF...\n")
                    download_cmd = (
                        f"aria2c -c -x 16 -s 16 -d {ref_dir} "
                        f"-o {genome_version}.gtf.gz {gtf_url}"
                    )
                    try:
                        subprocess.run(download_cmd, shell=True, check=True)
                        subprocess.run(f"gzip -d {gtf_gz}", shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        self.show_error_gui(f"GTF download failed: {str(e)}")
                        return
            else:
                self.update_output_gui("GTF file already exists.\n")

            # ====== tss file =========== (associated helper function is present after this method)
            tss_bed = os.path.join(ref_dir, f"{genome_version}_TSS.bed")
            if not os.path.exists(tss_bed):
                self.update_output_gui("Generating TSS BED file from GTF...\n")
                try:
                    self.gtf_to_tss_bed(
                        gtf_file=gtf_file,
                        output_bed=tss_bed,
                        genome_assembly=genome_version,
                        gene_biotype="protein_coding"
                    )
                except Exception as e:
                    self.show_error_gui(f"Failed to create TSS BED file: {str(e)}")
                    return

            # step 6: alignment and sorting
            self.update_output_gui("Checking Step 6: Alignment...\n")

            selected_samples = self.params["step1"].get("selected_samples")
            auto_sample_names = self.params["step1"].get("auto_sample_names", [])
            sample_file_map = {}

            if selected_samples:
                for r1_path in selected_samples:
                    base_r1 = os.path.basename(r1_path)

                    match = re.match(r"^(.*?)(_R?1(?:_001)?|_1)\.f(?:ast)?q\.gz$", base_r1)
                    if not match:
                        continue

                    sample_name = match.group(1)
                    r2_path = re.sub(r'(_R?)1(\d*)', r'\g<1>2\g<2>', r1_path)

                    if self.params["step2"].get("skip_trimming", False):
                        sample_file_map[sample_name] = {"r1": r1_path, "r2": r2_path}
                    else:
                        trimmed_r1 = os.path.join(trimmed_data, f"{sample_name}_val_1.fq.gz")
                        trimmed_r2 = os.path.join(trimmed_data, f"{sample_name}_val_2.fq.gz")

                        sample_file_map[sample_name] = {"r1": trimmed_r1, "r2": trimmed_r2}

                samples = list(sample_file_map.keys())

            else:
                samples = auto_sample_names
                for sample in samples:
                    if self.params["step2"].get("skip_trimming", False):
                        r1_candidates = glob.glob(os.path.join(raw_data, f"{sample}_*R1*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*1*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*R1*.fq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*1*.fq.gz"))

                        r2_candidates = glob.glob(os.path.join(raw_data, f"{sample}_*R2*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*2*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*R2*.fq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*2*.fq.gz"))

                        if not r1_candidates or not r2_candidates:
                            raise ValueError(f"Could not find paired files for sample {sample} in {raw_data}")

                        r1 = sorted(r1_candidates)[0]
                        r2 = sorted(r2_candidates)[0]

                    else:
                        r1 = os.path.join(trimmed_data, f"{sample}_val_1.fq.gz")
                        r2 = os.path.join(trimmed_data, f"{sample}_val_2.fq.gz")

                    sample_file_map[sample] = {"r1": r1, "r2": r2}

            # main alignment loop
            for sample in samples:
                coord_sorted_bam = os.path.join(bam_output, f"{sample}.sort.bam")
                name_sorted_bam = os.path.join(bam_output, f"{sample}.name.sorted.bam")

                if not (os.path.exists(coord_sorted_bam) and os.path.exists(coord_sorted_bam + ".bai")):
                    self.update_output_gui(f"Aligning {sample}...\n")

                    genome_version = self.params["step3"].get("genome_version", "").strip()

                    r1 = sample_file_map[sample]["r1"]
                    r2 = sample_file_map[sample]["r2"]

                    cmd = (f"bowtie2 -p {threads} --very-sensitive -X 2000 -x {bt2_base} "
                           f"-1 {r1} -2 {r2} "
                           f"| samtools view --threads {threads} -bS - > {bam_output}/{sample}.bam")
                    if not self.run_blocking_command(cmd):
                        return

                    if not os.path.exists(name_sorted_bam):
                        self.update_output_gui(f"Name sorting {sample}...\n")
                        cmd = f"samtools sort -n --threads {threads} -o {name_sorted_bam} {bam_output}/{sample}.bam && rm {bam_output}/{sample}.bam"
                        if not self.run_blocking_command(cmd):
                            return
                    else:
                        self.update_output_gui(f"Skipping name sorting for {sample} (Name-sorted BAM exists)\n")

                    self.update_output_gui(f"Coordinate sorting {sample}...\n")
                    cmd = f"samtools sort --threads {threads} -o {coord_sorted_bam} {name_sorted_bam} && samtools index -@ {threads} {coord_sorted_bam}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui(f"Skipping alignment for {sample} (BAM exists)\n")

            # flagstat
            flagstat_path = os.path.join(bam_output, "flagstat_results.txt")
            if os.path.exists(flagstat_path):
                os.remove(flagstat_path)

            bam_files = sorted(glob.glob(os.path.join(bam_output, "*.sort.bam")))

            with open(flagstat_path, "w") as stats_file:
                stats_file.write(
                    "Sample\tTotal Reads\tMapped (%)\tProperly Paired (%)\tSingletons (%)\tDiscordant (mapQ>=5)\n")

                for bam in bam_files:
                    name = os.path.basename(bam)
                    self.update_output_gui(f"Generating alignment stats for {name}...\n")

                    result = subprocess.run(
                        ["samtools", "flagstat", bam],
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    output = result.stdout

                    total_reads = re.search(r"(\d+) \+ \d+ in total", output)
                    mapped = re.search(r"(\d+) \+ \d+ mapped \(([\d\.]+)%", output)
                    properly_paired = re.search(r"(\d+) \+ \d+ properly paired \(([\d\.]+)%", output)
                    singletons = re.search(r"(\d+) \+ \d+ singletons \(([\d\.]+)%", output)
                    discordant = re.search(r"(\d+) \+ \d+ with mate mapped to a different chr \(mapQ>=5\)", output)

                    total = total_reads.group(1) if total_reads else "NA"
                    mapped_pct = mapped.group(2) if mapped else "NA"
                    proper_pct = properly_paired.group(2) if properly_paired else "NA"
                    single_pct = singletons.group(2) if singletons else "NA"
                    discordant_reads = discordant.group(1) if discordant else "NA"

                    stats_file.write(f"{name}\t{total}\t{mapped_pct}\t{proper_pct}\t{single_pct}\t{discordant_reads}\n")

            # step 7: coverage (associated helper function after this method)
            self.update_output_gui("Checking Step 7: Coverage...\n")
            for sample in samples:
                coverage_bw = os.path.join(normalized_coverage, f"{sample}.normalized.bw")
                if not os.path.exists(coverage_bw):
                    self.update_output_gui(f"Generating coverage for {sample}...\n")
                    cmd = f"bamCoverage -p {threads} -of bigwig --normalizeUsing=RPKM -b {bam_output}/{sample}.sort.bam -o {coverage_bw}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui(f"Skipping coverage for {sample} (file exists)\n")

            coverage_profile_dir = os.path.join(normalized_coverage, "coverage_profiles")
            os.makedirs(coverage_profile_dir, exist_ok=True)

            tss_prefix = os.path.join(coverage_profile_dir, "TSS")
            matrix_file = f"{tss_prefix}_matrix.gz"
            heatmap_file = f"{tss_prefix}_heatmap.pdf"

            if os.path.exists(matrix_file) and os.path.exists(heatmap_file):
                self.update_output_gui("TSS coverage profile and heatmap already exist. Skipping.\n")
            else:
                self.update_output_gui("Generating TSS coverage profiles...\n")
                try:
                    self.run_tss_matrix_and_heatmap(
                        tss_bed=tss_bed,
                        bw_dir=normalized_coverage,
                        out_prefix=tss_prefix
                    )
                except Exception as e:
                    self.show_error_gui(f"Failed to generate TSS matrix and heatmap: {str(e)}")
                    return

            # step 8: peak calling (genrich or macs3)

            # Genrich
            if self.params["step4"].get("peak_caller") == "Genrich":
                bam_dir = bam_output
                treated = self.params["step4"].get("treated_samples", [])
                control = self.params["step4"].get("control_samples", [])
                blacklist = ",".join(self.params["step4"].get("blacklist_files", []))
                exclude_chr_list = self.params["step4"].get("exclude_chr", [])
                if exclude_chr_list:
                    exclude_chr_arg = f"-e {','.join(exclude_chr_list)}"
                else:
                    exclude_chr_arg = ""
                expand_cut = self.params["step4"].get("expand_cut_sites", 100)
                max_q = self.params["step4"].get("max_q_value", 0.05)

                if blacklist:
                    blacklist_arg = f"-E {blacklist}"
                else:
                    blacklist_arg = ""

                treated_files = ",".join([f"{bam_dir}/{s}.name.sorted.bam" for s in treated])
                control_flag = f"-c {','.join([f'{bam_dir}/{s}.name.sorted.bam' for s in control])}" if control else ""

                if self.params["step4"].get("peak_type") == "merged":
                    output_peak = f"{peak_files}/{self.params['step4'].get('merged_output_name')}.genrich.peak"
                    output_bg = f"{peak_files}/{self.params['step4'].get('merged_output_name')}.genrich.peak.bg"
                    cmd_genrich = (f"Genrich -t {treated_files} {control_flag} -o {output_peak} -f {output_bg} "
                                   f"-v -j -d {expand_cut} -r -q {max_q} {exclude_chr_arg} {blacklist_arg}")
                    self.update_output_gui("Running Step 8: Genrich Peak Calling...\n")
                    if not self.run_blocking_command(cmd_genrich):
                        return
                    self.run_peak_annotation()

                else:
                    for sample in treated:
                        output_peak = f"{peak_files}/{sample}.genrich.peak"
                        output_bg = f"{peak_files}/{sample}.genrich.peak.bg"
                        cmd_genrich = (f"Genrich -t {bam_dir}/{sample}.name.sorted.bam {control_flag} "
                                       f"-o {output_peak} -f {output_bg} -v -j -d {expand_cut} -r "
                                       f"-q {max_q} {exclude_chr_arg} {blacklist_arg}")
                        self.update_output_gui(f"Running Genrich for {sample}...\n")
                        if not self.run_blocking_command(cmd_genrich):
                            return
                    self.run_peak_annotation()

            # MACS3
            else:
                bam_dir = bam_output
                treated = self.params["step4"].get("treated_samples", [])
                control = self.params["step4"].get("control_samples", [])
                blacklist_files = self.params["step4"].get("blacklist_files", [])
                exclude_chr = self.params["step4"].get("exclude_chr", [])
                max_q = self.params["step4"].get("max_q_value", 0.05)
                genome_size = self.params["step4"].get("genome", "hs")

                files_to_delete = set()
                control_filtered_map = {}
                macs3_commands = []

                def process_bam(sample, files_to_delete):
                    name_sorted_bam = os.path.join(bam_dir, f"{sample}.name.sorted.bam")
                    fixed_bam = os.path.join(bam_dir, f"{sample}.fixmate.bam")
                    coord_bam = os.path.join(bam_dir, f"{sample}.coord.bam")
                    dedup_bam = os.path.join(bam_dir, f"{sample}.dedup.bam")
                    filtered_bam = os.path.join(bam_dir, f"{sample}.filtered.dedup.bam")

                    cmd = (
                        f"samtools fixmate -m {name_sorted_bam} {fixed_bam} && "
                        f"samtools sort --threads {threads} -o {coord_bam} {fixed_bam} && "
                        f"samtools markdup -r --threads {threads} {coord_bam} {dedup_bam} && "
                        f"samtools index -@ {threads} {dedup_bam}"
                    )
                    files_to_delete.update([fixed_bam, coord_bam, dedup_bam])
                    return cmd, dedup_bam, filtered_bam

                #control processing
                for ctrl in control:
                    self.update_output_gui(f"Processing control sample {ctrl}...\n")
                    cmd, dedup, filtered = process_bam(ctrl, files_to_delete)
                    if not self.run_blocking_command(cmd): return
                    try:
                        self.filter_bam(dedup, filtered, blacklist_files, exclude_chr)
                        files_to_delete.add(filtered)
                        control_filtered_map[ctrl] = filtered
                    except Exception as e:
                        self.show_error_gui(f"BAM filtering failed for control sample {ctrl}: {e}")
                        return

                # process treated + construct macs3 command
                for s in treated:
                    self.update_output_gui(f"Processing test sample {s}...\n")
                    cmd, dedup, filtered = process_bam(s, files_to_delete)
                    if not self.run_blocking_command(cmd): return
                    try:
                        self.filter_bam(dedup, filtered, blacklist_files, exclude_chr)
                        files_to_delete.add(filtered)
                    except Exception as e:
                        self.show_error_gui(f"BAM filtering failed for sample {s}: {e}")
                        return

                    if control:
                        for ctrl in control:
                            control_bam = control_filtered_map[ctrl]
                            prefix = os.path.join(peak_files, f"{s}_{ctrl}.macs3.peak")
                            macs_cmd = (
                                f"macs3 callpeak -f BAMPE --call-summits -t {filtered} -c {control_bam} "
                                f"-g {genome_size} -n {prefix} -B -q {max_q}"
                            )
                            macs3_commands.append((f"{s} vs {ctrl}", macs_cmd))
                    else:
                        prefix = os.path.join(peak_files, f"{s}.macs3.peak")
                        macs_cmd = (
                            f"macs3 callpeak -f BAMPE --call-summits -t {filtered} "
                            f"-g {genome_size} -n {prefix} -B -q {max_q}"
                        )
                        macs3_commands.append((s, macs_cmd))

                self.update_output_gui("Running Step 8: MACS3 Peak Calling...\n")
                self.run_macs3_commands(
                    macs3_commands,
                    on_complete=lambda: self.cleanup_and_continue(files_to_delete)
                )

        except Exception as e:
            self.show_error_gui(f"Pipeline failed: {str(e)}")

    def run_macs3_commands(self, commands, on_complete=None):
        def worker():
            for label, cmd in commands:
                self.update_output_gui(f"\n Calling peaks for: {label}\n{cmd}\n")
                try:
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                            text=True)
                    for line in proc.stdout:
                        self.root.after(0, lambda l=line: self.update_output_gui(l))
                    exit_code = proc.wait()
                    if exit_code != 0:
                        self.root.after(0, lambda: self.show_error_gui(f"MACS3 failed for {label}"))
                        return
                except Exception as e:
                    self.root.after(0, lambda: self.show_error_gui(f"MACS3 crashed for {label}: {e}"))
                    return
            self.root.after(0, lambda: self.update_output_gui("MACS3 Peak Calling completed.\n"))
            if on_complete:
                self.root.after(0, on_complete)

        threading.Thread(target=worker, daemon=True).start()

    def cleanup_and_continue(self, files_to_delete):
        self.update_output_gui("Cleaning up temporary BAM files...\n")
        for f in files_to_delete:
            try:
                os.remove(f)
            except Exception as e:
                self.update_output_gui(f"Warning: Could not delete file {f}: {e}\n")
        self.run_peak_annotation()

    def run_peak_annotation(self):
        r_script_path = resource_filename('chromacs', 'annotate_peaks_5.R')
        if not os.path.exists(r_script_path):
            self.show_error_gui(f"R script 'annotate_peaks_5.R' not found at:\n{r_script_path}")
            return

        genome_version = self.params["step3"]["genome_version"]
        ref_dir = self.params["step3"]["ref_dir"]
        peak_files = os.path.join(self.params["step1"]["base_output_dir"], "peak_files")
        cmd = f"Rscript {r_script_path} {peak_files} {genome_version} {ref_dir}"
        self.update_output_gui(f"Running Peak Annotation with {r_script_path}...\n")

        def worker():
            try:
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        text=True)
                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))
                proc.wait()
                self.root.after(0, lambda: self.diffbind_btn.config(state=tk.NORMAL))
                timestamp = self.get_timestamp()
                self.update_output_gui(f"\n✅ Pipeline execution complete! [{timestamp}]\n")
                messagebox.showinfo("Pipeline Finished", "All steps have completed successfully!")
            except Exception as e:
                self.root.after(0, lambda: self.show_error_gui(f"Peak annotation failed: {e}"))

        threading.Thread(target=worker, daemon=True).start()

    # this is required for MACS3 filtering
    def filter_bam(self, input_bam, output_bam, blacklist_files, exclude_chr):
        import tempfile
        import shutil
        
        tmp_dir = tempfile.mkdtemp()
        step1 = os.path.join(tmp_dir, "step1.bam")

        try:
            # blacklist filtering
            if blacklist_files:
                blacklist_str = " ".join(blacklist_files)
                cmd1 = f"bedtools intersect -v -abam {input_bam} -b {blacklist_str} > {step1}"
            else:
                shutil.copy(input_bam, step1)
                cmd1 = None

            if cmd1:
                self.update_output_gui(f"Running blacklist filtering:\n{cmd1}\n")
                subprocess.run(cmd1, shell=True, check=True)

            # exclude chromosomes
            if exclude_chr:
                pat = "|".join(re.escape(c) for c in exclude_chr)
                cmd2 = (
                    f"samtools view -h {step1} | "
                    f"awk '$1 ~ /^@/ || $3 !~ /^({pat})$/' | "
                    f"samtools view -b - > {output_bam}"
                )
            else:
                shutil.copy(step1, output_bam)
                cmd2 = None

            if cmd2:
                self.update_output_gui(f"Running exclude chromosome filter:\n{cmd2}\n")
                subprocess.run(cmd2, shell=True, check=True)

        except subprocess.CalledProcessError as e:
            raise Exception(f"BAM filtering failed: {str(e)}")
        finally:
            shutil.rmtree(tmp_dir)

    # required for coverage matrix calculation and plotting of bigwigs
    def gtf_to_tss_bed(
            self,
            gtf_file: str,
            output_bed: str,
            genome_assembly: str,
            gene_biotype: str = "protein_coding"
    ):
        mito_chromosomes = {
            "GRCh38": {"chrM", "MT", "M"},
            "GRCm39": {"chrM", "MT", "M"},
            "mRatBN7.2": {"chrM", "MT", "M"},
            "ARS-UCD1.3": {"chrM", "MT", "M"},  # Bos taurus mitochondrion is usually "MT" or "chrM"
            "Sscrofa11.1": {"chrM", "MT", "M"},
            "GRCg7b": {"chrM", "MT", "M"},  # Gallus gallus (chicken) mito often "chrM" or "MT"
            "Pan_tro_3.0": {"chrM", "MT", "M"},  # Chimpanzee mitochondrion likely similar to human
            "ROS_Cfam_1.0": {"chrM", "MT", "M"},  # Dog mitochondrion
            "ARS1": {"chrM", "MT", "M"},  # Capra hircus (goat) mitochondrion, check assembly docs if needed
            "CVASU_BBG_1.0": {"chrM", "MT", "M"},  # Another goat assembly, same as above
            "OryCun2.0": {"chrM", "MT", "M"},  # Rabbit mitochondrion
            "gorGor4": {"chrM", "MT", "M"},  # Gorilla mitochondrion
            "Mmul_10": {"chrM", "MT", "M"},  # Macaca mulatta (rhesus macaque)

            "GRCz11": {"chrM", "MT", "M"},
            "UCB_Xtro_10.0": {"chrM", "MT", "M"},  # Xenopus tropicalis
            "Ssal_v3.1": {"chrM", "MT", "M"},  # Salmo salar (Atlantic salmon)

            "BDGP6": {"chrM", "mitochondrion", "MT"},
            "WBcel235": {"MtDNA", "chrM"},
        }

        # Keyword-based blacklist
        non_primary_keywords = {"random", "alt", "hap", "fix", "patch", "scaffold", "gl", "un", "ki"}

        def is_primary_contig(chrom: str) -> bool:
            chrom = chrom.lower()
            return not any(k in chrom for k in non_primary_keywords)

        exclude_chroms = mito_chromosomes.get(genome_assembly, set())

        with (
                gzip.open(gtf_file, "rt") if gtf_file.endswith(".gz") else open(gtf_file, "r")
        ) as infile, open(output_bed, "w") as outfile:

            seen_genes = set()

            for line in infile:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 9 or cols[2].lower() != "gene":
                    continue

                chrom = cols[0]

                if chrom in exclude_chroms:
                    continue

                if not is_primary_contig(chrom):
                    continue

                attrs = {}
                for pair in cols[8].strip().split(";"):
                    pair = pair.strip()
                    if not pair:
                        continue
                    if "\"" in pair:
                        key, val = pair.split(" ", 1)
                        val = val.strip("\"")
                    else:
                        key, val = pair.split("=") if "=" in pair else (pair, "")
                    attrs[key.strip()] = val.strip()

                biotype = attrs.get("gene_biotype") or attrs.get("gene_type") or attrs.get("biotype")
                if not biotype or biotype.lower() != gene_biotype.lower():
                    continue

                strand = cols[6]
                start = int(cols[3])
                end = int(cols[4])
                tss = start if strand == "+" else end

                bed_start = tss - 1
                bed_end = bed_start + 1

                gene_id = attrs.get("gene_id") or attrs.get("ID")
                if not gene_id or gene_id in seen_genes:
                    continue
                seen_genes.add(gene_id)

                outfile.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_id}\t.\t{strand}\n")

    def run_tss_matrix_and_heatmap(self, tss_bed: str, bw_dir: str, out_prefix: str):

        bw_files = glob.glob(os.path.join(bw_dir, "*.normalized.bw"))
        if not bw_files:
            raise RuntimeError(f"No BigWig files found in {bw_dir}")

        matrix_file = f"{out_prefix}_matrix.gz"
        cmd_mat = [
                      "computeMatrix", "reference-point",
                      "--referencePoint", "center",
                      "-b", "3000", "-a", "3000",
                      "-R", tss_bed,
                      "-S"
                  ] + bw_files + [
                      "--skipZeros",
                      "-p", "max",
                      "-o", matrix_file
                  ]
        subprocess.run(cmd_mat, check=True)

        heatmap_file = f"{out_prefix}_heatmap.pdf"
        cmd_heat = [
            "plotHeatmap",
            "-m", matrix_file,
            "-out", heatmap_file,
            "--colorMap", "Reds",
            "--heatmapHeight", "15",
            "--heatmapWidth", "20",
            "--xAxisLabel", "Distance from TSS (bp)",
            "--refPointLabel", "TSS",
            "--averageTypeSummaryPlot", "mean",
            "--plotType", "se",
            "--legendLocation", "best"
        ]
        subprocess.run(cmd_heat, check=True)

        return matrix_file, heatmap_file

    # =================== diffbind analysis logic =====================================================================

    def launch_diffbind_config(self):
        self.diffbind_window = tk.Toplevel(self.root)
        self.diffbind_window.title("DiffBind Configuration")
        self.diffbind_window.geometry("1000x700")

        self.all_samples = self.params["step2"].get("sample_names", [])

        tk.Label(self.diffbind_window, text="1. Select Peak Files:", font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, padx=10, pady=5)
        self.peak_listbox = tk.Listbox(self.diffbind_window, selectmode=tk.MULTIPLE, exportselection=False, height=10,width=80)
        self.peak_listbox.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.diffbind_window, text="Browse Peaks",command=self.populate_peak_list, bg='gray').grid(row=0, column=2)

        self.config_frame = tk.Frame(self.diffbind_window)
        self.config_frame.grid(row=1, column=0, columnspan=3, padx=10, pady=10, sticky='nw')

        instruction_1 = (
            "Note: You can skip this Control section, if there is none\n"
            "Please keep in mind that these are not biological controls, rather technical controls\n"
            "Do NOT put baseline condition as control here\n\n"
            "IMPORTANT: “Input the baseline and test conditions accordingly\n"
            "Think baseline as your reference, and test as the comparison group\n"
            "These labels can represent any two biological states — e.g.:\n"
            "Baseline = Lung tissue, Test = Kidney tissue\n"
            "Baseline = Wild-type, Test = Mutant\n"
            "Baseline = Day 0, Test = Day 5\n"
        )

        tk.Label(
            self.diffbind_window, text=instruction_1,
            wraplength=650, justify=tk.LEFT, anchor="w",
            font=(self.roboto_font, 9)
        ).grid(row=1, column=2, padx=10, sticky="w")

        tk.Label(self.diffbind_window, text="FDR Threshold (0–1):", font=(self.roboto_font, 10)).grid(row=2, column=0,padx=10,sticky='w')
        self.fdr_threshold = tk.DoubleVar(value=0.05)
        tk.Spinbox(self.diffbind_window, from_=0.0, to=1.0, increment=0.01, textvariable=self.fdr_threshold,width=5).grid(row=2, column=1, sticky='w')


        tk.Label(self.diffbind_window, text="Number of Threads (default = 4):", font=(self.roboto_font, 10)).grid(row=4, column=0, sticky='w', padx=10)
        self.diffbind_threads = tk.StringVar(value="4")  # default is 4
        tk.Entry(self.diffbind_window, textvariable=self.diffbind_threads, width=5).grid(row=4, column=1, sticky='w')

        instruction_1 = (
            "Note: DiffBind is very memory intensive tool, so adjust your thread here accordingly\n"
            "      If run fails with error in the core setup, adjust (decrease thread count) and rerun DiffBind\n"
        )
        tk.Label(
            self.diffbind_window, text=instruction_1,
            wraplength=650, justify=tk.LEFT, anchor="w",
            font=(self.roboto_font, 9)
        ).grid(row=6, column=1, padx=10, sticky="w")

        tk.Button(self.diffbind_window, text="Run DiffBind", command=self.validate_and_run_diffbind, bg="yellow green").grid(row=8, column=1, pady=10)

        # metadata scores
        self.conditions = {}
        self.replicates = {}
        self.controls = {}

        self.peak_listbox.bind('<<ListboxSelect>>', lambda e: self.build_param_rows())

        tk.Button(
            self.diffbind_window,
            text="Annotate Results",
            command=self.launch_diffbind_annotation,
            bg="cyan"
        ).grid(row=9, column=1, pady=10)

    def populate_peak_list(self):
        self.peak_listbox.delete(0, tk.END)
        peak_dir = os.path.join(self.params["step1"]["base_output_dir"], "peak_files")
        peak_files = glob.glob(os.path.join(peak_dir, "*.narrowPeak")) \
                     + glob.glob(os.path.join(peak_dir, "*.genrich.peak"))
        for pf in sorted(peak_files):
            name = os.path.basename(pf)
            self.peak_listbox.insert(tk.END, name)
            self.conditions[name] = tk.StringVar(value="test") # default
            self.replicates[name] = tk.IntVar(value=1)
            self.controls[name] = tk.StringVar(value="")

    def build_param_rows(self):
        for child in self.config_frame.winfo_children():
            child.destroy()

        headers = ['Peak File', 'Condition', 'Replicate', 'Control']
        for col, h in enumerate(headers):
            tk.Label(self.config_frame, text=h, font=(self.roboto_font, 9, 'bold')).grid(row=0, column=col, padx=5)

        for i, idx in enumerate(self.peak_listbox.curselection(), start=1):
            peak = self.peak_listbox.get(idx)
            if peak not in self.conditions:
                self.conditions[peak] = tk.StringVar(value="test")
            if peak not in self.replicates:
                self.replicates[peak] = tk.IntVar(value=1)
            if peak not in self.controls:
                self.controls[peak] = tk.StringVar(value="")

            tk.Label(self.config_frame, text=peak).grid(row=i, column=0, sticky='w')

            ttk.Combobox(self.config_frame, textvariable=self.conditions[peak], values=["test", "baseline"], state="readonly").grid(row=i, column=1)

            tk.Spinbox(self.config_frame, from_=1, to=10,
                       textvariable=self.replicates[peak], width=5).grid(row=i, column=2)

            ttk.Combobox(self.config_frame, textvariable=self.controls[peak],
                         values=self.all_samples, state="readonly").grid(row=i, column=3)

    def validate_and_run_diffbind(self):
        if not self.peak_listbox.curselection():
            return self.show_error_gui("Please select at least one peak file.")

        has_any_control = False
        missing_controls = []

        for idx in self.peak_listbox.curselection():
            peak = self.peak_listbox.get(idx)
            if not self.conditions[peak].get():
                return self.show_error_gui(f"Condition missing for {peak}")
            if not self.replicates[peak].get():
                return self.show_error_gui(f"Replicate missing for {peak}")
            if self.controls[peak].get():
                has_any_control = True
            else:
                missing_controls.append(peak)

        if has_any_control and missing_controls:
            peak_list = "\n".join(missing_controls)
            proceed = messagebox.askyesno(
                "Mixed Control Usage",
                f"Some peak files are missing control inputs:\n\n{peak_list}\n\n"
                "Do you want to continue without controls for those?"
            )
            if not proceed:
                return

        fdr = self.fdr_threshold.get()
        if not (0.0 <= fdr <= 1.0):
            return self.show_error_gui("FDR threshold must be between 0 and 1.")

        self.run_diffbind_analysis()

    def generate_metadata(self):
        rows = []
        base = self.params["step1"]["base_output_dir"]
        has_control = False

        for idx in self.peak_listbox.curselection():
            peak_filename = self.peak_listbox.get(idx)
            sample_id_with_suffix = os.path.basename(peak_filename).split('.')[0]

            if 'macs3' in peak_filename.lower():
                sample_id = sample_id_with_suffix.split('_')[0]
                peak_caller = 'narrow'
            else:
                sample_id = sample_id_with_suffix
                peak_caller = 'bed'

            condition = self.conditions[peak_filename].get()
            replicate = self.replicates[peak_filename].get()
            control_id = self.controls[peak_filename].get()
            if control_id: has_control = True

            bam_reads = os.path.join(base, "bam_output", f"{sample_id}.sort.bam")
            peak_path = os.path.join(base, "peak_files", peak_filename)
            bam_control = os.path.join(base, "bam_output", f"{control_id}.sort.bam") if control_id else ""

            row = {
                "SampleID": sample_id,
                "bamReads": bam_reads,
                "Condition": condition,
                "Replicate": replicate,
                "Peaks": peak_path,
                "PeakCaller": peak_caller
            }

            if control_id:
                row["ControlID"] = control_id
                row["bamControl"] = bam_control

            rows.append(row)

        df = pd.DataFrame(rows)

        if has_control:
            cols = ["SampleID", "bamReads", "Condition", "Replicate",
                    "Peaks", "PeakCaller", "ControlID", "bamControl"]
        else:
            cols = ["SampleID", "bamReads", "Condition", "Replicate",
                    "Peaks", "PeakCaller"]

        return df[cols]

    def run_diffbind_analysis(self):

        base_dir = self.params["step1"]["base_output_dir"]
        out_dir = os.path.join(base_dir, "diffbind_results")
        os.makedirs(out_dir, exist_ok=True)

        df = self.generate_metadata()
        out_csv = os.path.join(out_dir, "diffbind_samplesheet.csv")
        df.to_csv(out_csv, index=False)

        r_script_path = resource_filename('chromacs', 'diffbind_3.R')
        fdr = str(self.fdr_threshold.get())
        threads = self.diffbind_threads.get().strip()
        threads = threads if threads.isdigit() and int(threads) > 0 else "4"  # fallback to 4 if invalid

        cmd = ["Rscript", r_script_path, out_csv, out_dir, fdr, threads]

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )

                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))

                exit_code = proc.wait()

                if exit_code == 0:
                    def on_success():
                        timestamp = self.get_timestamp()
                        self.update_output_gui(f"DiffBind analysis completed successfully! [{timestamp}]\n")
                        messagebox.showinfo(
                            "DiffBind Completed",
                            f"Results are in:\n{out_dir}"
                        )

                    self.root.after(0, on_success)

                else:
                    def on_error():
                        msg = f" DiffBind failed with exit code {exit_code}\n"
                        self.update_output_gui(msg)
                        self.show_error_gui(msg)

                    self.root.after(0, on_error)

            except Exception as e:
                def on_exception():
                    msg = f" Unexpected error running DiffBind:\n{e}\n"
                    self.update_output_gui(msg)
                    self.show_error_gui(msg)

                self.root.after(0, on_exception)

        threading.Thread(target=worker, daemon=True).start()

    def launch_diffbind_annotation(self):
        base_dir = self.params["step1"]["base_output_dir"]
        diffbind_dir = os.path.join(base_dir, "diffbind_results")
        ref_dir = self.params["step3"]["ref_dir"]
        assembly = self.params["step3"]["genome_version"]

        r_script = resource_filename('chromacs', 'annotate_diffbind.R')

        cmd = [
            "Rscript", r_script,
            diffbind_dir,
            assembly,
            ref_dir
        ]

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))
                exit_code = proc.wait()

                if exit_code != 0:
                    msg = f" Annotate DiffBind script failed with exit code {exit_code}\n"
                    self.root.after(0, lambda: self.show_error_gui(msg))
                    self.root.after(0, lambda: self.update_output_gui(msg))
                    return

                try:
                    self._inject_fimo_peakids(
                        annotated_dir=os.path.join(diffbind_dir, "Annotated_DiffBind"),
                        fimo_bed_path=os.path.join(diffbind_dir, "consensus_peaks_for_fimo.bed"),
                        output_column_name="peak_id"
                    )

                except Exception as e:
                    err_msg = f"Failed injecting FIMO IDs: {e}\n"
                    self.root.after(0, lambda: self.update_output_gui(err_msg))
                    self.root.after(0, lambda: self.show_error_gui(err_msg))
                    return

                self.root.after(0, lambda: messagebox.showinfo(
                    "Success", "DiffBind annotation complete!"
                ))
            except Exception as e:
                self.root.after(0, lambda: self.show_error_gui(str(e)))

        threading.Thread(target=worker, daemon=True).start()

    def _inject_fimo_peakids(self, annotated_dir, fimo_bed_path, output_column_name="peak_id"):
        if not os.path.exists(annotated_dir):
            raise FileNotFoundError(f"Annotated directory not found: {annotated_dir}")

        if not os.path.exists(fimo_bed_path):
            self.root.after(0, lambda: self.update_output_gui(f"FIMO bed not found: {fimo_bed_path}\n"))
            return

        fimo_df = pd.read_csv(
            fimo_bed_path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start", "end", "peakID"],
            dtype={"chrom": str, "start": int, "end": int, "peakID": str}
        )

        coord_map = defaultdict(list)
        for _, row in fimo_df.iterrows():
            chrom = row.chrom
            start = row.start
            end = row.end
            pid = row.peakID
            for ds in (-1, 0, 1):
                for de in (-1, 0, 1):
                    key = (chrom, start + ds, end + de)
                    coord_map[key].append(pid)

        for annotated_file in os.listdir(annotated_dir):
            if not annotated_file.endswith("_annotated.tsv"):
                continue
            file_path = os.path.join(annotated_dir, annotated_file)
            try:
                ann_df = pd.read_csv(file_path, sep="\t", dtype=str)
            except Exception as e:
                self.root.after(0, lambda: self.update_output_gui(f"Failed reading {annotated_file}: {e}\n"))
                continue

            cols_lower = {c.lower(): c for c in ann_df.columns}
            chr_col = next((cols_lower[n] for n in ("seqnames", "chr", "chrom", "chromosome") if n in cols_lower), None)
            start_col = next((cols_lower[n] for n in ("start", "chromstart") if n in cols_lower), None)
            end_col = next((cols_lower[n] for n in ("end", "chromend") if n in cols_lower), None)

            if chr_col is None or start_col is None or end_col is None:
                self.root.after(0, lambda: self.update_output_gui(
                    f"Skipping {annotated_file}: missing chr/start/end columns\n"))
                continue

            fimo_ids = []
            matched = 0
            unmatched = 0
            for _, row in ann_df.iterrows():
                chrom = str(row[chr_col])
                try:
                    start = int(row[start_col])
                    end = int(row[end_col])
                except (ValueError, TypeError):
                    fimo_ids.append(None)
                    unmatched += 1
                    continue

                key = (chrom, start, end)
                peakids = coord_map.get(key)
                if not peakids:
                    peakids = coord_map.get((chrom, start, end), None)

                if peakids:
                    unique_ids = sorted(set(peakids))
                    fimo_ids.append(";".join(unique_ids))
                    matched += 1
                else:
                    fimo_ids.append(None)
                    unmatched += 1

            colname = output_column_name
            if colname in ann_df.columns:
                ann_df[colname] = fimo_ids
            else:
                insert_pos = 3 if ann_df.shape[1] >= 3 else ann_df.shape[1]
                ann_df.insert(insert_pos, colname, fimo_ids)

            summary = f"{annotated_file}: matched {matched}, unmatched {unmatched}\n"
            self.root.after(0, lambda s=summary: self.update_output_gui(s))

            ann_df = self._merge_with_diffbind_noisq_metrics(annotated_df=ann_df, annotated_dir=annotated_dir,
                                                             annotated_file=annotated_file, peakid_column=colname)

            ann_df.to_csv(file_path, sep="\t", index=False)

    def _merge_with_diffbind_noisq_metrics(self, annotated_df, annotated_dir, annotated_file, peakid_column="peak_id"):
        ann_df = annotated_df
        parent_dir = os.path.dirname(annotated_dir)
        results_df = None
        results_file = None

        diffbind_results = os.path.join(parent_dir, "diffbind_results.csv")
        if os.path.exists(diffbind_results):
            results_file = diffbind_results
            results_df = pd.read_csv(results_file)
        else:
            noisq_results = os.path.join(os.path.dirname(parent_dir), "noisq_results", "noisq_results.xlsx")
            if os.path.exists(noisq_results):
                results_file = noisq_results
                results_df = pd.read_excel(results_file)

        if results_df is None:
            self.root.after(0, lambda: self.update_output_gui(
                f"No diffbind/noisq results found to merge for {annotated_file}\n"))
            return ann_df

        peakid_col = None
        for candidate in ("PeakID", "peak_id", "peakID"):
            if candidate in results_df.columns:
                peakid_col = candidate
                break
        if peakid_col is None:
            self.root.after(0, lambda: self.update_output_gui(
                f"Results file {os.path.basename(results_file)} lacks a peak ID column; skipping merge for {annotated_file}\n"))
            return ann_df

        other_cols = [c for c in results_df.columns if c != peakid_col]
        last_six = other_cols[-6:] if len(other_cols) >= 6 else other_cols
        subset = results_df[[peakid_col] + last_six].drop_duplicates()

        ann_df["merge_key"] = ann_df[peakid_column].fillna("").str.split(";").str[0]

        ann_df["merge_key"] = ann_df[peakid_column].fillna("").astype(str).str.split(";").str[0]
        subset[peakid_col] = subset[peakid_col].astype(str)

        subset = subset.rename(columns={peakid_col: "results_peak_id"})

        merged = ann_df.merge(
            subset,
            left_on="merge_key",
            right_on="results_peak_id",
            how="left"
        )

        merged.drop(columns=["merge_key", "results_peak_id"], inplace=True, errors="ignore")

        #reproting (optional)
        matched = 0
        if last_six:
            matched = merged[last_six[0]].notna().sum()
        #self.root.after(0, lambda: self.update_output_gui(f"{annotated_file}: merged {matched}/{len(ann_df)} rows from {os.path.basename(results_file)}\n"))

        return merged

    # ======================  noiseq analysis logic =====================================================================

    def launch_noisq_config(self):
        self.noisq_window = tk.Toplevel(self.root)
        self.noisq_window.title("NOISeq Configuration")
        self.noisq_window.geometry("1000x600")

        tk.Label(self.noisq_window, text="Select Peak Files (1 per sample):", font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, padx=10, pady=5)

        self.noisq_peak_listbox = tk.Listbox(self.noisq_window, selectmode=tk.MULTIPLE, exportselection=False,height=10, width=80)
        self.noisq_peak_listbox.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.noisq_window, text="Browse Peaks",command=self.populate_noisq_peak_list, bg='gray').grid(row=0, column=2)


        self.noisq_meta_frame = tk.Frame(self.noisq_window)
        self.noisq_meta_frame.grid(row=1, column=0, columnspan=2, padx=(10, 40), pady=10, sticky="nw")

        self.noisq_conditions = {}

        self.noisq_peak_listbox.bind('<<ListboxSelect>>', lambda e: self.build_noisq_param_rows())

        instruction_1 = (
            "IMPORTANT: “Prove baseline and test condition accordingly\n"
            "Think baseline as your reference, and test as the comparison group\n"
            "These labels can represent any two biological states — e.g.:\n"
            "Baseline = Lung tissue, Test = Kidney tissue\n"
            "Baseline = Wild-type, Test = Mutant\n"
            "Baseline = Day 0, Test = Day 5\n"
        )

        tk.Label(
            self.noisq_window, text=instruction_1,
            wraplength=350, justify=tk.LEFT, anchor="n",
            font=(self.roboto_font, 9)
        ).grid(row=1, column=2, padx=(0, 10), sticky="nw")

        tk.Label(self.noisq_window, text="Confidence Threshold (NOISeq q, 0–1):",font=(self.roboto_font, 10)).grid(row=2, column=0, padx=10, sticky='w')
        self.noisq_qvalue = tk.DoubleVar(value=0.9)
        tk.Spinbox(self.noisq_window, from_=0.0, to=1.0, increment=0.01,textvariable=self.noisq_qvalue, width=5).grid(row=2, column=1, sticky='w')

        tk.Button(self.noisq_window, text="Run NOISeq",command=self.run_noisq_gui, bg="light green").grid(row=3, column=1, pady=10)

        tk.Button(self.noisq_window, text="Annotate NOISeq Results", command=self.launch_noisq_annotation,bg="cyan").grid(row=4, column=1, pady=10)

    def populate_noisq_peak_list(self):
        self.noisq_peak_listbox.delete(0, tk.END)
        peak_dir = os.path.join(self.params["step1"]["base_output_dir"], "peak_files")
        peak_files = glob.glob(os.path.join(peak_dir, "*.narrowPeak")) \
                     + glob.glob(os.path.join(peak_dir, "*.genrich.peak"))
        for pf in sorted(peak_files):
            name = os.path.basename(pf)
            self.noisq_peak_listbox.insert(tk.END, name)
            self.noisq_conditions[name] = tk.StringVar(value="test")

    def build_noisq_param_rows(self):
        for child in self.noisq_meta_frame.winfo_children():
            child.destroy()

        tk.Label(self.noisq_meta_frame, text="Peak File").grid(row=0, column=0, padx=5)
        tk.Label(self.noisq_meta_frame, text="Condition").grid(row=0, column=1, padx=5)

        for i, idx in enumerate(self.noisq_peak_listbox.curselection(), start=1):
            peak = self.noisq_peak_listbox.get(idx)
            tk.Label(self.noisq_meta_frame, text=peak).grid(row=i, column=0, sticky='w')
            ttk.Combobox(self.noisq_meta_frame, textvariable=self.noisq_conditions[peak],
                         values=["test", "baseline"], state="readonly").grid(row=i, column=1)

    def generate_noisq_metadata(self):
        rows = []
        base = self.params["step1"]["base_output_dir"]

        for idx in self.noisq_peak_listbox.curselection():
            peak_filename = self.noisq_peak_listbox.get(idx)
            sample_id_with_suffix = os.path.basename(peak_filename).split('.')[0]

            if 'macs3' in peak_filename.lower():
                sample_id = sample_id_with_suffix.split('_')[0]
            else:
                sample_id = sample_id_with_suffix

            peak_path = os.path.join(base, "peak_files", peak_filename)
            bam_reads = os.path.join(base, "bam_output", f"{sample_id}.sort.bam")

            condition = self.noisq_conditions[peak_filename].get()

            rows.append({
                "SampleID": sample_id,
                "Condition": condition,
                "Peaks": peak_path,
                "bamReads": bam_reads,
            })

        return pd.DataFrame(rows)

    def run_noisq_gui(self):
        selected = self.noisq_peak_listbox.curselection()
        if not selected:
            self.show_error_gui("Please select at least one peak file for NOISeq.")
            return

        base_dir = self.params["step1"]["base_output_dir"]
        peak_dir = os.path.join(base_dir, "peak_files")
        out_dir = os.path.join(base_dir, "noisq_results")
        os.makedirs(out_dir, exist_ok=True)

        metadata_df = self.generate_noisq_metadata()
        meta_csv = os.path.join(out_dir, "metadata.csv")
        metadata_df.to_csv(meta_csv, index=False)
        threads = self.params.get("threads", 8)

        # consensus peak
        try:
            all_peaks = [os.path.join(peak_dir, self.noisq_peak_listbox.get(idx)) for idx in selected]
            merged_bed = os.path.join(out_dir, "consensus_peaks.bed")

            cat_cmd = f"cat {' '.join(all_peaks)} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > {merged_bed}"
            result = subprocess.run(cat_cmd, shell=True, check=True, capture_output=True, text=True)
            self.update_output_gui(result.stdout)

            #saf file
            saf_file = os.path.join(out_dir, "consensus_peaks.saf")
            bed4_file = os.path.join(out_dir, "consensus_peaks_for_fimo.bed")
            with open(merged_bed, 'r') as infile, \
                    open(saf_file, 'w') as saf_out, \
                    open(bed4_file, 'w') as bed4_out:
                saf_out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
                for i, line in enumerate(infile, 1):
                    parts = line.strip().split('\t')
                    peak_id = f"peak{i}"

                    saf_out.write(f"{peak_id}\t{parts[0]}\t{parts[1]}\t{parts[2]}\t.\n")
                    bed4_out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{peak_id}\n")

            #featurecounts
            bam_files = metadata_df["bamReads"].tolist()
            count_matrix = os.path.join(out_dir, "peak_counts.txt")
            feature_cmd = (
                f"featureCounts -a {saf_file} -F SAF -T {threads} -p -B -C "
                f"-o {count_matrix} {' '.join(bam_files)}"
            )
            result = subprocess.run(feature_cmd, shell=True, check=True, capture_output=True, text=True)
            self.update_output_gui(result.stdout)

            counts_df = pd.read_csv(count_matrix, sep='\t', comment='#', skiprows=1)
            counts_df = counts_df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
            counts_df = counts_df.rename(columns={"Geneid": "peak_id"})

            bam_to_sample = dict(zip(metadata_df["bamReads"], metadata_df["SampleID"]))
            new_columns = []
            for col in counts_df.columns:
                if col == "peak_id":
                    new_columns.append(col)
                else:
                    new_columns.append(bam_to_sample.get(col, col))
            counts_df.columns = new_columns

            peak_matrix_clean = os.path.join(out_dir, "peak_matrix.csv")
            counts_df.to_csv(peak_matrix_clean, index=False)

        except subprocess.CalledProcessError as e:
            self.show_error_gui(f"Command failed: {e.cmd}\nError: {e.stderr}")
            return

        r_script = resource_filename('chromacs', 'noisq_atac.R')
        output_xlsx = os.path.join(out_dir, "noisq_results.xlsx")
        counts_csv = os.path.join(out_dir, "peak_matrix.csv")
        qval = str(self.noisq_qvalue.get())
        cmd = ["Rscript", r_script, counts_csv, meta_csv, output_xlsx, qval]

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                for line in iter(proc.stdout.readline, ''):
                    self.update_output_gui(line)
                proc.wait()
                if proc.returncode == 0:
                    timestamp = self.get_timestamp()
                    self.update_output_gui(f"NOISeq analysis completed! [{timestamp}]\n")
                    messagebox.showinfo("Success", f"Results saved to:\n{output_xlsx}")
                else:
                    self.show_error_gui(f"NOISeq failed with exit code {proc.returncode}")
            except Exception as e:
                self.show_error_gui(f"Error running NOISeq: {str(e)}")

        threading.Thread(target=worker, daemon=True).start()

    def launch_noisq_annotation(self):
        base_dir = self.params["step1"]["base_output_dir"]
        noisq_results = os.path.join(base_dir, "noisq_results", "noisq_results.xlsx")
        peak_counts = os.path.join(base_dir, "noisq_results", "peak_counts.txt")
        assembly = self.params["step3"]["genome_version"]
        ref_dir = self.params["step3"]["ref_dir"]

        r_script = resource_filename('chromacs', 'annotate_noisq.R')

        cmd = [
            "Rscript", r_script,
            noisq_results,
            peak_counts,
            assembly,
            ref_dir
        ]

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))
                exit_code = proc.wait()

                if exit_code != 0:
                    msg = f" Annotate NOISeq script failed with exit code {exit_code}\n"
                    self.root.after(0, lambda: self.show_error_gui(msg))
                    self.root.after(0, lambda: self.update_output_gui(msg))
                    return

                # Post-process FIMO IDs injection
                try:
                    noisq_dir = os.path.join(base_dir, "noisq_results")
                    self._inject_fimo_peakids(
                        annotated_dir=os.path.join(noisq_dir, "Annotated_NOISeq"),
                        fimo_bed_path=os.path.join(noisq_dir, "consensus_peaks_for_fimo.bed"),
                        output_column_name="peak_id"
                    )


                except Exception as e:
                    err_msg = f"Failed injecting FIMO IDs: {e}\n"
                    self.root.after(0, lambda: self.update_output_gui(err_msg))
                    self.root.after(0, lambda: self.show_error_gui(err_msg))
                    return

                self.root.after(0, lambda: messagebox.showinfo(
                    "Success", "NOISeq annotation complete!"
                ))


            except Exception as e:
                self.root.after(0, lambda: self.show_error_gui(str(e)))

        threading.Thread(target=worker, daemon=True).start()


    # ==================== motif enrichment logic =====================================================================
    def launch_motif_module(self):
        self.motif_window = tk.Toplevel(self.root)
        self.motif_window.title("Motif Enrichment")
        self.motif_window.geometry("1000x600")

        self.motif_window.grid_columnconfigure(1, weight=1)
        for i in range(8):
            self.motif_window.grid_rowconfigure(i, weight=0)

        title_label = tk.Label(
            self.motif_window,
            text="Motif Enrichment Configuration",
            font=(self.roboto_font, 11, "bold")
        )
        title_label.grid(row=0, column=0, columnspan=3, pady=10, sticky="n")

        file_selectors = [
            ("Differential Peaks File (CSV/XLSX/TSV)", "motif_diff_file", 1),
            ("Merged Regions BED File (consensus_peaks_for_fimo)", "motif_bed_file", 2),
            ("Genome FASTA File", "motif_fasta_file", 3),
            ("Motif MEME File", "motif_meme_file", 4)
        ]

        for label_text, attr, row in file_selectors:
            tk.Label(
                self.motif_window,
                text=label_text,
                anchor="w"
            ).grid(row=row, column=0, padx=10, pady=5, sticky="w")

            entry = tk.Entry(self.motif_window, width=80)
            entry.grid(row=row, column=1, padx=10, pady=5, sticky="ew")
            setattr(self, f"{attr}_entry", entry)

            tk.Button(
                self.motif_window,
                text="Browse",
                command=lambda a=attr: self.browse_file(a)
            ).grid(row=row, column=2, padx=10, pady=5, sticky="e")

        genome_version = self.params["step3"].get("genome_version")
        ref_dir = self.params["step3"].get("ref_dir")
        if genome_version and ref_dir:
            fasta_path = os.path.join(ref_dir, f"{genome_version}.fa")
            if os.path.exists(fasta_path):
                self.motif_fasta_file_entry.delete(0, tk.END)
                self.motif_fasta_file_entry.insert(0, fasta_path)

        tk.Label(
            self.motif_window,
            text="Output Directory:",
            anchor="w"
        ).grid(row=5, column=0, padx=10, pady=5, sticky="w")

        base_dir = self.params["step1"]["base_output_dir"]
        motif_out_dir = os.path.join(base_dir, "motif_results")
        self.output_dir_label = tk.Label(
            self.motif_window,
            text=motif_out_dir,
            relief="sunken",
            anchor="w",
            bg="white"
        )
        self.output_dir_label.grid(row=5, column=1, padx=10, pady=5, sticky="ew", columnspan=2)

        run_button = tk.Button(
            self.motif_window,
            text="Run Motif Enrichment Pipeline",
            command=self.run_motif_enrichment_pipeline,
            bg="light green",
            padx=20,
            pady=10
        )
        run_button.grid(row=7, column=0, columnspan=3, pady=20, sticky="n")

        self.motif_window.grid_rowconfigure(8, weight=1)

    def run_motif_enrichment_pipeline(self):
        base_dir = self.params["step1"]["base_output_dir"]
        motif_out_dir = os.path.join(base_dir, "motif_results")
        os.makedirs(motif_out_dir, exist_ok=True)

        diff_file = self.motif_diff_file_entry.get()
        bed_file = self.motif_bed_file_entry.get()
        fasta_file = self.motif_fasta_file_entry.get()
        meme_file = self.motif_meme_file_entry.get()
        output_prefix = os.path.join(motif_out_dir, "motif")

        def worker():
            try:
                import pandas as pd
                from scipy.stats import fisher_exact
                from statsmodels.stats.multitest import multipletests
                import numpy as np
                import os

                self._update_output("Starting motif enrichment pipeline...\n")

                # Check if output enrichment files exist
                up_enrichment_file = f"{output_prefix}_up_enrichment.tsv"
                down_enrichment_file = f"{output_prefix}_down_enrichment.tsv"

                if os.path.exists(up_enrichment_file) and os.path.exists(down_enrichment_file):
                    self._update_output(
                        "Enrichment result files detected. Skipping analysis and going to plotting...\n")

                    self._update_output("Loading enrichment files...\n")

                    try:
                        up_df = pd.read_csv(up_enrichment_file, sep='\t')
                        down_df = pd.read_csv(down_enrichment_file, sep='\t')
                        print("Loaded files successfully!")
                    except Exception as e:
                        self._update_output(f"Failed to load enrichment files: {str(e)}\n")
                        raise


                else:

                    fasta_output = f"{output_prefix}_peak_seqs.fa"
                    self._update_output("Generating FASTA sequences from BED file...\n")
                    self.run_blocking_command(f"bedtools getfasta -fi '{fasta_file}' -bed '{bed_file}' -fo '{fasta_output}' -name")

                    fimo_output = f"{output_prefix}_motifs.fimo.tsv"
                    self._update_output("Running FIMO motif scanning...\n")
                    self.run_blocking_command(f"fimo --text --skip-matched-sequence '{meme_file}' '{fasta_output}' > '{fimo_output}'", show_output=False)

                    self._update_output("Running motif enrichment analysis...\n")


                    if diff_file.endswith('.csv'):
                        diff_df = pd.read_csv(diff_file)
                    elif diff_file.endswith('.xlsx'):
                        diff_df = pd.read_excel(diff_file)
                    else:
                        diff_df = pd.read_csv(diff_file, sep='\t')

                    fold_col_candidates = ["Fold", "M"]
                    peakid_col_candidates = ["PeakID", "peak_id"]
                    fold_col = next((col for col in fold_col_candidates if col in diff_df.columns), None)
                    peakid_col = next((col for col in peakid_col_candidates if col in diff_df.columns), None)
                    if not fold_col or not peakid_col:
                        raise ValueError("Required columns not found in differential peaks file.")
                    up_peaks = set(diff_df.loc[diff_df[fold_col] > 0, peakid_col].astype(str))
                    down_peaks = set(diff_df.loc[diff_df[fold_col] < 0, peakid_col].astype(str))

                    motif_map = self.parse_meme_motif_names(meme_file)

                    fimo_df = pd.read_csv(fimo_output, sep='\t')
                    fimo_df.columns = [col.lstrip('#') for col in fimo_df.columns]  # remove '#' from header
                    fimo_df = fimo_df[['pattern name', 'sequence name', 'p-value']]
                    fimo_df.columns = ['motif_id', 'sequence_name', 'p-value']
                    fimo_df['sequence_name'] = fimo_df['sequence_name'].str.replace(r'::.*', '', regex=True)

                    bed_df = pd.read_csv(bed_file, sep='\t', header=None)
                    if bed_df.shape[1] < 4:
                        raise ValueError("BED file must have at least 4 columns with peak IDs in the 4th column.")
                    columns = ['chr', 'start', 'end', 'peak_id', 'score', 'strand'][:bed_df.shape[1]]
                    bed_df.columns = columns
                    bed_df['peak_id'] = bed_df['peak_id'].astype(str).str.strip('"').str.strip("'")
                    all_peaks = set(bed_df['peak_id'])

                    diff_df[diff_df.columns[0]] = diff_df[diff_df.columns[0]].astype(str).str.strip('"').str.strip("'")
                    fimo_df['sequence_name'] = fimo_df['sequence_name'].astype(str).str.strip('"').str.strip("'")
                    bed_df['peak_id'] = bed_df['peak_id'].astype(str).str.strip('"').str.strip("'")

                    for group_name, group_peaks in [('up', up_peaks), ('down', down_peaks)]:
                        self._update_output(f"Calculating enrichment for {group_name}-regulated peaks...\n")

                        group_hits = fimo_df[fimo_df['sequence_name'].isin(group_peaks)]
                        background_hits = fimo_df[fimo_df['sequence_name'].isin(all_peaks)]

                        results = []
                        for motif_id in fimo_df['motif_id'].unique():
                            group_peaks_with_motif = set(group_hits[group_hits['motif_id'] == motif_id]['sequence_name'])
                            background_peaks_with_motif = set(background_hits[background_hits['motif_id'] == motif_id]['sequence_name'])

                            A = len(group_peaks_with_motif)  # group peaks with motif
                            B = len(group_peaks) - A  # group peaks without motif

                            bg_peaks_excluding_group = all_peaks - group_peaks

                            #motif hits in background only (exclude group peaks)
                            C = len(background_peaks_with_motif.intersection(bg_peaks_excluding_group))
                            D = len(bg_peaks_excluding_group) - C  # background peaks without motif

                            if min(A, B, C, D) < 0:
                                continue

                            odds, pvalue = fisher_exact([[A, C], [B, D]], alternative='greater')
                            percent_group = A / (A + B) if (A + B) > 0 else 0
                            percent_background = C / (C + D) if (C + D) > 0 else 0
                            log_odds = np.log2(odds) if odds > 0 else 0

                            motif_name = motif_map.get(motif_id, motif_id)
                            results.append({
                                'Motif': motif_id,
                                'Motif_Name': motif_name,
                                'Hits_in_diff': A,
                                'Unhit_in_diff': B,
                                'Hits_in_background': C,
                                'Unhit_in_background': D,
                                'Percent_in_diff': percent_group,
                                'Percent_in_background': percent_background,
                                'Log2_Odds_Ratio': log_odds,
                                'PValue': pvalue
                            })

                        result_df = pd.DataFrame(results)

                        if not result_df.empty:
                            result_df['QValue'] = multipletests(result_df['PValue'], method='fdr_bh')[1]
                            result_df = result_df.sort_values('PValue')
                        else:
                            self._update_output(f"No enriched motifs found for {group_name} peaks.\n")

                        output_file = f"{output_prefix}_{group_name}_enrichment.tsv"
                        result_df.to_csv(output_file, sep='\t', index=False)
                        self._update_output(f"Saved {group_name} enrichment results to {output_file}\n")

                # plotting
                self._update_output("Generating visualizations...\n")

                up_df = pd.read_csv(f"{output_prefix}_up_enrichment.tsv", sep='\t')
                down_df = pd.read_csv(f"{output_prefix}_down_enrichment.tsv", sep='\t')

                self.plot_top_motifs_bar(up_df, "Up", os.path.join(motif_out_dir, "motif_top_up_bar.png"))
                self.plot_top_motifs_bar(down_df, "Down", os.path.join(motif_out_dir, "motif_top_down_bar.png"))
                self.plot_volcano(up_df, "Up", os.path.join(motif_out_dir, "motif_up_volcano.png"))
                self.plot_volcano(down_df, "Down", os.path.join(motif_out_dir, "motif_down_volcano.png"))
                self.plot_enrichment_heatmap(up_df, down_df, os.path.join(motif_out_dir, "motif_enrichment_heatmap.png"))

                self._update_output("Creating PDF report...\n")
                plot_files = [
                    os.path.join(motif_out_dir, "motif_top_up_bar.png"),
                    os.path.join(motif_out_dir, "motif_top_down_bar.png"),
                    os.path.join(motif_out_dir, "motif_up_volcano.png"),
                    os.path.join(motif_out_dir, "motif_down_volcano.png"),
                    os.path.join(motif_out_dir, "motif_enrichment_heatmap.png")
                ]
                self.create_pdf_report(plot_files, os.path.join(motif_out_dir, "motif_enrichment_report.pdf"))

                plot_output_dir = os.path.join(motif_out_dir, "motif_gene_analysis")
                os.makedirs(plot_output_dir, exist_ok=True)

                self._update_output("Starting motif annotation...\n")
                self.annotate_motif_results()

                plot_output_dir = os.path.join(motif_out_dir, "motif_gene_analysis")
                os.makedirs(plot_output_dir, exist_ok=True)

                self._update_output("Loading json files...\n")
                up_df = self.load_motif_data(os.path.join(motif_out_dir, "motif_up_annotated.json"))
                down_df = self.load_motif_data(os.path.join(motif_out_dir, "motif_down_annotated.json"))

                self._update_output("Plotting for upregulated...\n")
                self.plot_top_motifs_by_genes(up_df, "Upregulated", os.path.join(plot_output_dir, "up_motifs_by_genes.png"))
                self.plot_motif_gene_network(up_df, os.path.join(plot_output_dir, "up_motif_gene_network.png"))
                self.plot_gene_expression_heatmap(up_df, "Upregulated", os.path.join(plot_output_dir, "up_gene_heatmap.png"))

                self._update_output("Plotting for downregulated...\n")
                self.plot_top_motifs_by_genes(down_df, "Downregulated", os.path.join(plot_output_dir, "down_motifs_by_genes.png"))
                self.plot_motif_gene_network(down_df, os.path.join(plot_output_dir, "down_motif_gene_network.png"))
                self.plot_gene_expression_heatmap(down_df, "Downregulated",
                                             os.path.join(plot_output_dir, "down_gene_heatmap.png"))

                self._update_output("Plotting for combined...\n")
                combined_df = pd.concat([up_df, down_df])
                self.plot_top_motifs_by_genes(combined_df, "All Significant",
                                         os.path.join(plot_output_dir, "combined_motifs_by_genes.png"))

                plot_files.extend([
                    os.path.join(plot_output_dir, "up_motifs_by_genes.png"),
                    os.path.join(plot_output_dir, "up_motif_gene_network.png"),
                    os.path.join(plot_output_dir, "up_gene_heatmap.png"),
                    os.path.join(plot_output_dir, "down_motifs_by_genes.png"),
                    os.path.join(plot_output_dir, "down_motif_gene_network.png"),
                    os.path.join(plot_output_dir, "down_gene_heatmap.png"),
                    os.path.join(plot_output_dir, "combined_motifs_by_genes.png")
                ])

                self._update_output("Plotting of motif:gene complete\n")
                self.show_info_gui("Motif enrichment pipeline completed successfully!")

            except Exception as e:
                self.show_error_gui(f"Error in motif enrichment pipeline:\n{str(e)}")

        threading.Thread(target=worker, daemon=True).start()

    def browse_file(self, attr):
        if attr == "motif_diff_file":
            filetypes = [
                ("Differential Peaks", "*.csv *.xlsx"),
                ("All files", "*")
            ]
        elif attr == "motif_bed_file":
            filetypes = [("BED files", "*.bed"), ("All files", "*")]
        elif attr == "motif_fasta_file":
            filetypes = [("FASTA files", "*.fa *.fasta"), ("All files", "*")]
        elif attr == "motif_meme_file":
            filetypes = [("MEME files", "*.meme *.txt"), ("All files", "*")]
        else:
            filetypes = [("All files", "*")]

        file_path = filedialog.askopenfilename(parent=self.motif_window, filetypes=filetypes)

        if file_path:
            entry = getattr(self, f"{attr}_entry")
            entry.delete(0, tk.END)
            entry.insert(0, file_path)


    def parse_meme_motif_names(self, meme_file):
        motif_map = {}
        with open(meme_file) as f:
            for line in f:
                if line.startswith('MOTIF'):
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        motif_id = parts[1]
                        motif_name = parts[2]
                        # uniqueness as motif name may be multiple times
                        motif_label = f"{motif_name} ({motif_id})"
                        motif_map[motif_id] = motif_label
                    elif len(parts) == 2:
                        motif_id = parts[1]
                        motif_map[motif_id] = motif_id
        return motif_map
    
#-------------------------------- motif enrichment plots ---------------------------------------------
    def plot_top_motifs_bar(self, df, group_name, out_file, top_n=20):
        self._update_output(f"Generating {group_name} bar plot...\n")
        df = df.copy()
        df['-log10(QValue)'] = -np.log10(df['QValue'] + 1e-10)
        df = df.sort_values('QValue').head(top_n)
        plt.figure(figsize=(10, 6))
        sns.barplot(x='-log10(QValue)', y='Motif_Name', data=df, hue='Motif_Name', palette='viridis', legend=False)
        plt.title(f'Top {top_n} Enriched Motifs ({group_name})')
        plt.xlabel('-log10(Q-value)')
        plt.tight_layout()
        plt.savefig(out_file)
        plt.close()

    def plot_volcano(self, df, group_name, out_file):
        self._update_output(f"Generating {group_name} volcano plot...\n")
        df = df.copy()
        df['-log10(PValue)'] = -np.log10(df['PValue'] + 1e-10)
        plt.figure(figsize=(8, 6))
        sns.scatterplot(data=df, x='Log2_Odds_Ratio', y='-log10(PValue)',
                        hue='QValue', palette='coolwarm', legend=False)
        plt.axhline(-np.log10(0.05), color='gray', linestyle='--')
        plt.title(f'Motif Enrichment Volcano Plot ({group_name})')
        plt.xlabel('Log2 Odds Ratio')
        plt.ylabel('-log10(P-value)')
        plt.tight_layout()
        plt.savefig(out_file)
        plt.close()

    def plot_enrichment_heatmap(self, up_df, down_df, out_file):
        self._update_output(f"Generating {up_df} vs {down_df} heatmap...\n")
        top_up = up_df.nsmallest(20, 'QValue')
        top_down = down_df.nsmallest(20, 'QValue')
        top_motifs = set(top_up['Motif_Name']).union(set(top_down['Motif_Name']))

        up_df = up_df[up_df['Motif_Name'].isin(top_motifs)].set_index('Motif_Name')
        down_df = down_df[down_df['Motif_Name'].isin(top_motifs)].set_index('Motif_Name')

        heatmap_df = pd.DataFrame({
            'Up': up_df['Log2_Odds_Ratio'],
            'Down': down_df['Log2_Odds_Ratio']
        }).fillna(0)

        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_df, annot=True, cmap='vlag', center=0)
        plt.title('Motif Enrichment Heatmap (Log2 OR)')
        plt.ylabel('Motif')
        plt.xlabel('Group')
        plt.tight_layout()
        plt.savefig(out_file)
        plt.close()

    def create_pdf_report(self, image_paths, pdf_path):
        self._update_output(f"Generating pdf from images...\n")
        time.sleep(4)
        
        pdf = FPDF()
        for img_path in image_paths:
            cover = Image.open(img_path)
            width, height = cover.size
            pdf_w, pdf_h = 210, 297  # A4 in mm
            img_w = pdf_w
            img_h = height * (pdf_w / width)

            pdf.add_page()
            pdf.image(img_path, x=0, y=0, w=img_w, h=img_h)
        pdf.output(pdf_path, "F")
#-------------------------------------------------------------------------------------------------------

    def annotate_motif_results(self):
        base_dir = self.params["step1"]["base_output_dir"]
        motif_out_dir = os.path.join(base_dir, "motif_results")
        diff_file = self.motif_diff_file_entry.get()

        try:
            fimo_file = os.path.join(motif_out_dir, "motif_motifs.fimo.tsv")
            fimo_df = pd.read_csv(fimo_file, sep='\t')
            fimo_df.columns = [col.lstrip('#').strip() for col in fimo_df.columns]

            fimo_df[['peak_id', 'coord_key']] = fimo_df['sequence name'].str.split('::', expand=True)
            fimo_df['coord_key'] = fimo_df['coord_key'].str.replace('_', ':')

            analysis_type = "diffbind" if "diffbind" in diff_file.lower() else "noisq"

            def load_annotations(analysis_type):
                if analysis_type == "diffbind":
                    gain_path = os.path.join(base_dir, "diffbind_results", "Annotated_DiffBind",
                                             "gain_sites_annotated.tsv")
                    loss_path = os.path.join(base_dir, "diffbind_results", "Annotated_DiffBind",
                                             "loss_sites_annotated.tsv")
                else:
                    gain_path = os.path.join(base_dir, "noisq_results", "Annotated_NOISeq", "gain_sites_annotated.tsv")
                    loss_path = os.path.join(base_dir, "noisq_results", "Annotated_NOISeq", "loss_sites_annotated.tsv")

                try:
                    gain_df = pd.read_csv(gain_path, sep='\t')
                    loss_df = pd.read_csv(loss_path, sep='\t')

                    def create_key(row, adjust_start=False):
                        chr_col = next((col for col in ['seqnames', 'chr'] if col in row.index), None)
                        if not chr_col:
                            raise ValueError("Could not find chromosome column in annotation file")

                        start = row['start'] - 1 if adjust_start else row['start']
                        return f"{row[chr_col]}:{start}-{row['end']}"

                    adjust_start = (analysis_type == "diffbind")
                    gain_df['coord_key'] = gain_df.apply(create_key, axis=1, adjust_start=adjust_start)
                    loss_df['coord_key'] = loss_df.apply(create_key, axis=1, adjust_start=adjust_start)

                    return gain_df, loss_df

                except Exception as e:
                    self._update_output(f"Error loading annotation files: {str(e)}\n")
                    raise

            gain_df, loss_df = load_annotations(analysis_type)

            def process_motif_file(motif_file, annotation_df, output_suffix):
                try:
                    motif_df = pd.read_csv(motif_file, sep='\t')
                    sig_motifs = motif_df[motif_df['QValue'] < 0.05]

                    if sig_motifs.empty:
                        self._update_output(f"No significant motifs found for {output_suffix}\n")
                        return

                    motif_fimo = fimo_df[fimo_df['pattern name'].isin(sig_motifs['Motif'])]

                    if motif_fimo.empty:
                        self._update_output(f"No FIMO hits found for significant {output_suffix} motifs\n")
                        return

                    merged = pd.merge(
                        motif_fimo,
                        annotation_df,
                        left_on='coord_key',
                        right_on='coord_key',
                        how='left',
                        indicator=True
                    )

                    matched_fraction = (merged['_merge'] == 'both').mean()
                    self._update_output(f"True annotation match rate: {matched_fraction:.1%}\n")

                    records = []
                    for motif_name, grp in merged.groupby('pattern name', dropna=False):
                        gene_cols = ['SYMBOL', 'geneId', 'GENENAME', 'annotation', 'distanceToTSS']
                        available_gene_cols = [c for c in gene_cols if c in grp.columns]
                        if available_gene_cols:
                            genes_df = grp[available_gene_cols].drop_duplicates()
                            genes = genes_df.to_dict('records')
                        else:
                            genes = []

                        peaks = sorted(grp['peak_id'].dropna().unique().tolist())
                        peak_count = len(peaks)

                        records.append({
                            'Motif': motif_name,
                            'Annotations': {
                                'genes': genes,
                                'peak_count': peak_count,
                                'peaks': peaks
                            }
                        })

                    motif_annotations = pd.DataFrame(records)

                    final_df = pd.merge(
                        sig_motifs,
                        motif_annotations,
                        left_on='Motif',
                        right_on='Motif',
                        how='left'
                    )

                    if 'Motif_Name' not in final_df.columns:
                        final_df['Motif_Name'] = final_df['Motif']

                    output_file = os.path.join(motif_out_dir, f"motif_{output_suffix}_annotated.json")
                    final_df.to_json(output_file, orient='records', indent=2)
                    self._update_output(f"Saved annotated {output_suffix} motifs to {output_file}\n")

                except Exception as e:
                    self._update_output(f"Error processing {output_suffix} motifs: {str(e)}\n")
                    raise

            process_motif_file(
                os.path.join(motif_out_dir, "motif_up_enrichment.tsv"),
                gain_df,
                "up"
            )
            process_motif_file(
                os.path.join(motif_out_dir, "motif_down_enrichment.tsv"),
                loss_df,
                "down"
            )

            self._update_output("Motif annotation completed successfully!\n")

        except Exception as e:
            self.show_error_gui(f"Error in motif annotation: {str(e)}")

    def load_motif_data(self, json_path):
        with open(json_path) as f:
            data = json.load(f)

        processed = []
        for motif in data:
            genes = [g for g in motif['Annotations']['genes'] if g['SYMBOL'] is not None]
            gene_count = len(genes)
            gene_symbols = [g['SYMBOL'] for g in genes if g['SYMBOL']]

            processed.append({
                'Motif': motif['Motif'],
                'Motif_Name': motif['Motif_Name'],
                'QValue': motif['QValue'],
                'Log2_Odds_Ratio': motif['Log2_Odds_Ratio'],
                'Gene_Count': gene_count,
                'Genes': gene_symbols,
                'Peak_Count': motif['Annotations']['peak_count'],
                'Top_Genes': ', '.join(list(set(gene_symbols))[:5])
            })

        return pd.DataFrame(processed)

#----------------------------------------------------motif annotation plots----------------------------------------
    def plot_top_motifs_by_genes(self, df, title, output_path, top_n=20):
        top_motifs = df.sort_values('Gene_Count', ascending=False).head(top_n)

        plt.figure(figsize=(12, 8))
        sns.barplot(data=top_motifs, x='Gene_Count', y='Motif_Name', hue='Motif_Name',
                    palette='viridis', dodge=False, legend=False)

        plt.title(f'Top {top_n} {title} Motifs by Number of Target Genes')
        plt.xlabel('Number of Target Genes')
        plt.ylabel('Motif')
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()

    def plot_motif_gene_network(self, df, output_path, top_motifs=10, top_genes=20):
        top_motifs_df = df.sort_values('Gene_Count', ascending=False).head(top_motifs)
        all_genes = [g for sublist in top_motifs_df['Genes'] for g in sublist]
        gene_counts = Counter(all_genes)
        top_genes = [gene for gene, count in gene_counts.most_common(top_genes)]

        edges = []
        for _, row in top_motifs_df.iterrows():
            for gene in row['Genes']:
                if gene in top_genes:
                    edges.append((row['Motif_Name'], gene))

        G = nx.Graph()
        G.add_edges_from(edges)

        plt.figure(figsize=(15, 10))
        pos = nx.spring_layout(G, k=0.5)

        nx.draw_networkx_nodes(G, pos, nodelist=top_motifs_df['Motif_Name'].tolist(),
                               node_color='lightblue', node_size=2000, label='Motifs')
        nx.draw_networkx_nodes(G, pos, nodelist=top_genes,
                               node_color='lightgreen', node_size=1000, label='Genes')
        nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
        nx.draw_networkx_labels(G, pos, font_size=10)

        plt.title(f'Network of Top {top_motifs} Motifs and Their Target Genes')
        plt.legend()
        plt.axis('off')
        plt.savefig(output_path)
        plt.close()

    def plot_gene_expression_heatmap(self, df, title, output_path):
        all_genes = [g for sublist in df['Genes'] for g in sublist]
        gene_counts = pd.Series(all_genes).value_counts().head(20)

        plt.figure(figsize=(10, 8))
        sns.barplot(x=gene_counts.values, y=gene_counts.index, palette='rocket')
        plt.title(f'Most Frequent Target Genes in {title} Motifs')
        plt.xlabel('Number of Motifs Binding This Gene')
        plt.ylabel('Gene Symbol')
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

    def run_blocking_command(self, command, show_output=True):
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        output_lines = []

        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                output_lines.append(output)
                if show_output:
                    self._update_output(output)

        returncode = process.wait()

        if returncode != 0:
            error_summary = "".join(output_lines[-10:]) if not show_output else ""
            self.show_error_gui(f"Command failed: {command}\n{error_summary}")
            return False

        return True

    def update_output_gui(self, text):
        self.root.after(0, self._update_output, text)

    def _update_output(self, text):
        self.output_text.config(state=tk.NORMAL)
        self.output_text.insert(tk.END, text)
        self.output_text.see(tk.END)
        self.output_text.config(state=tk.DISABLED)

    def show_error_gui(self, message):
        self.root.after(0, lambda: messagebox.showerror("Error", message))

    def show_info_gui(self, message):
        self.root.after(0, lambda: messagebox.showinfo("Info", message))


def main():
    root = tk.Tk()

    icon_path = resource_filename("chromacs", "assets/ChromAcS.png")
    if os.path.exists(icon_path):
        try:
            img = tk.PhotoImage(file=icon_path)
            root.tk.call('wm', 'iconphoto', root._w, img)
        except Exception as e:
            print(f"Could not set app icon: {e}")

    app = ATACSeqPipeline(root)
    root.mainloop()


if __name__ == "__main__":
    main()

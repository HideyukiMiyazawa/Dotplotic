import tkinter as tk
from tkinter import filedialog, messagebox
import subprocess
import threading

def choose_blast_file():
    filename = filedialog.askopenfilename(filetypes=[("BLAST file", "*")])
    blast_entry.delete(0, tk.END)
    blast_entry.insert(0, filename)

def choose_annotation_file():
    filename = filedialog.askopenfilename(filetypes=[("Annotation file", "*")])
    annotation_entry.delete(0, tk.END)
    annotation_entry.insert(0, filename)

def choose_output_file():
    filename = filedialog.asksaveasfilename(filetypes=[("SVG files", "Dotplotic.svg")])
    output_entry.delete(0, tk.END)
    output_entry.insert(0, filename)

def run_perl_script():
    # Disable
    run_button.config(state="disabled")
    status_label.config(text="Running...", justify=['right'])
    root.update()

    # Processing
    threading.Thread(target=run_perl_script_worker).start()
    
def run_perl_script_worker():
    blast_file = blast_entry.get()
    annotation_file = annotation_entry.get()
    output_file = output_entry.get()
    query = q_entry.get()
    sbjct = s_entry.get()
    outfmt = outfmt_entry.get()
    IdnRng = IdnRng_entry.get()
    Ann = annotation_entry.get()
    Colset = colset_var.get()
    Bgcolset = bg_colset_var.get()

    if not blast_file or not output_file:
        messagebox.showwarning("Error")
        return

    if layout_var.get():
        Layout = "on"
    else:
        Layout = "off"

    command = ["perl", "Dotplotic", 
        "--blast", blast_file, 
        "--out", output_file, 
        "--query", query, 
        "--subject", sbjct,
        "--outfmt", outfmt,
        "--identity_range", IdnRng,
        "--annotation", Ann,
        "--color_set", Colset,
        "--bg_color_set", Bgcolset,
        "--auto_layout", Layout,
    ]

    if not al_var.get():
        command.append("--align_direction")

    if click_var.get():
        command.append("--click")

    if light_var.get():
        command.append("--light")

    try:
        subprocess.run(command, check=True)
        messagebox.showinfo("Success", "Dotplotic finished successfully.")
    except subprocess.CalledProcessError as e:
        messagebox.showerror("Failed", f"Dotplotic failed\n{e}")
    finally:
        # Update bottons and labels
        status_label.config(text="")
        run_button.config(state="normal")

# building of GUI
root = tk.Tk()
root.title("Dotplotic GUI")

# input BLAST file
tk.Label(root, text="BLAST file:").grid(row=0, column=0, sticky="e")
blast_entry = tk.Entry(root, width=25)
blast_entry.insert(0, "test_data/Blast.tsv")
blast_entry.grid(row=0, column=1)
tk.Button(root, text="select", command=choose_blast_file).grid(row=0, column=2, sticky="w")

# query
tk.Label(root, text="Query:").grid(row=1, column=0, sticky="e")
q_entry = tk.Entry(root, width=25)
q_entry.grid(row=1, column=1)

# subject
tk.Label(root, text="Subject:").grid(row=2, column=0, sticky="e")
s_entry = tk.Entry(root, width=25)
s_entry.grid(row=2, column=1)

# BLAST outfmt
tk.Label(root, text="BLAST outfmt:").grid(row=3, column=0, sticky="e")
outfmt_entry = tk.Entry(root, width=25)
outfmt_entry.insert(0, "6 std")
outfmt_entry.grid(row=3, column=1)

# Identity range (def: 60-100)
tk.Label(root, text="Identity range:").grid(row=4, column=0, sticky="e")
IdnRng_entry = tk.Entry(root, width=25)
IdnRng_entry.insert(0, "60-100")
IdnRng_entry.grid(row=4, column=1)

# Annotation file
tk.Label(root, text="Annotation file:").grid(row=5, column=0, sticky="e")
annotation_entry = tk.Entry(root, width=25)
annotation_entry.insert(0, "")
annotation_entry.grid(row=5, column=1)
tk.Button(root, text="select", command=choose_annotation_file).grid(row=5, column=2, sticky="w")

# align_direction
al_var = tk.BooleanVar(value=True)
al_bool = tk.Checkbutton(root, text="Auto align direction", variable=al_var)
al_bool.grid(row=6, column=1, sticky='w', padx=10, pady=5)

# click
click_var = tk.BooleanVar(value=False)
click_bool = tk.Checkbutton(root, text="Show info on click", variable=click_var)
click_bool.grid(row=7, column=1, sticky='w', padx=10, pady=5)

# light version
light_var = tk.BooleanVar(value=False)
light_bool = tk.Checkbutton(root, text="Light-version", variable=light_var)
light_bool.grid(row=8, column=1, sticky='w', padx=10, pady=5)

# color_set
colset_var = tk.StringVar(value="1")
tk.Label(root, text="Alingment color set:").grid(row=9, column=0, sticky="e")
tk.Radiobutton(root, text="Red, Yellow, Light Green", variable=colset_var, value="1").grid(row=9, column=1, sticky="w")
tk.Radiobutton(root, text="Blue, Cyan, Yellow-Orange", variable=colset_var, value="2").grid(row=9, column=2, sticky="w")

# bg_color_set
bg_colset_var = tk.StringVar(value="1")
tk.Label(root, text="Background color:").grid(row=10, column=0, sticky="e")
tk.Radiobutton(root, text="Black", variable=bg_colset_var, value="1").grid(row=10, column=1, sticky="w")
tk.Radiobutton(root, text="White", variable=bg_colset_var, value="2").grid(row=10, column=2, sticky="w")

# auto layout
layout_var = tk.BooleanVar(value=True)
layout_bool = tk.Checkbutton(root, text="Auto layout", variable=layout_var)
layout_bool.grid(row=11, column=1, sticky='w', padx=10, pady=5)

# output SVG file name
tk.Label(root, text="").grid(row=12, column=0, pady=10)
tk.Label(root, text="output SVG file:").grid(row=13, column=0, sticky="e")
output_entry = tk.Entry(root, width=25)
output_entry.insert(0, "Dotplotic.svg")
output_entry.grid(row=13, column=1)
tk.Button(root, text="select", command=choose_output_file).grid(row=13, column=2, sticky="w")

# Run & Quit buttons
tk.Label(root, text="").grid(row=14, column=0, pady=10)
button_frame = tk.Frame(root)
button_frame.grid(row=15, column=0, columnspan=3, pady=10)
run_button = tk.Button(button_frame, text="Run", command=run_perl_script)
run_button.pack(side="left", padx=10)
quit_button = tk.Button(button_frame, text="Quit", command=root.destroy)
quit_button.pack(side="left", padx=10)

status_label = tk.Label(root, text="")
status_label.grid(row=16, pady=10)

#
# Main
#

if __name__ == "__main__":
    root.mainloop()
# Tutorials for Using Computing Systems and Simulations  
*Hall Group – Graduate Student Onboarding*

This tutorial provides a **practical introduction** to computing systems and simulation workflows used in the Hall Group.  
It is intended as a **reference guide**, not something to memorize. New students should skim it once and return to relevant sections as needed.

---

## How to Use This Tutorial

- You are **not expected to understand everything on the first read**
- Use this document when you:
  - Forget a command
  - Need to compile or debug a code
  - Submit jobs to HPC
  - Visualize simulation output
- Italicized text indicates placeholders that depend on your file names  
- Commands appear in `code blocks`

---

## 1. Basic Computing Tools (Concepts, Not Software)

Regardless of operating system or software choice, you should be comfortable with:

- Editing **plain-text files**
- Transferring files between **local and remote systems**
- Using a **terminal (SSH)** to run programs
- Compiling and debugging simulation codes

Suggested tools include Notepad++, VS Code, WinSCP, PuTTY, or native terminal tools on macOS/Linux.

---

## 2. File Types You Will Encounter

### Text Files
Human-readable files containing code, inputs, or outputs.

Examples: `.txt .out .err .cpp .f90 .py .readme`


### Executables
Compiled programs that can be run but not read as text.

Examples: `.exe .bin`


### Compressed Files
Collections of files bundled together.

Examples: `.zip .tar`


> Most simulation codes begin as text files and are compiled into executables.

---

## 3. Editing Files

### Local Editing
Use a text editor with:
- Line numbers
- Syntax highlighting

These features are essential for debugging.

### Remote Editing
- Small edits: edit directly via file-transfer tools
- Larger edits: download → edit locally → upload

⚠️ Avoid opening very large output files in GUI tools — use terminal commands instead.

---

## 4. Accessing Remote Systems (SSH)

You will use SSH to:
- Navigate directories
- Compile and run codes
- Submit HPC jobs

After logging in, you will see a **command prompt** indicating the system is ready for input.

---

## 5. Essential Terminal Commands

| Task | Command |
|----|----|
| Show current directory | `pwd` |
| List directory contents | `ls` |
| Change directory | `cd path/` |
| Go up one level | `cd ../` |
| Run executable | `./executable_name` |
| Stop a command | `Ctrl + C` |
| Log out | `logout` |

> Commands are **case-sensitive**.

You can scroll through previous commands using the **up arrow**, and autocomplete file names using **Tab**.

---

## 6. Running Executables

Executables must be run from their directory:

```bash
./executable_name
```

With required arguments:

```bash
./executable_name arg1 arg2 arg3
```

If unsure which arguments are required, run the executable without arguments and read the error message.

### Permission Errors

If you see “permission denied”:

```bash
chmod +x executable_name
```

## 7. Transferring Files

Files can be transferred using:

- Drag-and-drop via GUI tools such as WinSCP

- Terminal commands

After running simulations, refresh directories to see newly generated files.

⚠️ Large files should be inspected using terminal tools or analysis scripts.

## 8. Group Workstations

We currently have 1 CPU and 2 GPUs workstations for the exclusive use of members of the Hall group. If you are new to our group, contact the OIT Help Desk to get access to our workstations.

## 9. High-Performance Computing (HPC)

HPC access: `login.hpc.ncsu.edu`

Your working directory: `/gpfs_common/share/username/`

Important Rules

❌ Do not run simulations on login nodes
✅ Use login nodes only for editing, compiling, and job submission

Jobs run on compute nodes via the scheduling system.

### Checking Job Status - bjobs

Useful options:
```
bjobs -u all -q queue_name
bjobs -u all -q queue_name | less
```

### Submitting Jobs to HPC

Jobs are submitted using a script file. Example Submission Script:

```
#!/bin/bash
#BSUB -q single_chassis
#BSUB -n 2
#BSUB -W 120
#BSUB -J job_name
#BSUB -o sim%J.out
#BSUB -e sim%J.err

./executable_name input_file
```

Submit with: `bsub < script_name`

Jobs may wait in the queue before starting.

### Extracting Compressed Files

To extract .tar or .zip files: `tar -zxf filename`

Extracted files will appear in the current directory.

### Editing Files with Vi

Vi is a terminal-based text editor.

Common Commands:

| Action | Command |
|------|--------|
| Open file | `vi file_name` |
| Enter insert mode | `i` |
| Save and quit | `:wq` |
| Quit | `:q` |
| Quit without saving | `:q!` |

Press `Esc` to exit insert mode.

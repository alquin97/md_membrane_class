# Protocol

Before you start, clone this repo anywhere you want in your local machine.

```
git clone https://github.com/alquin97/md_membrane_class
```

## Simple membrane system

### Just POPC

This part of the protocol will be done under `just_popc/`:

```
cd just_popc
```

#### Membrane creation with [PACKMOL-Memgen](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00269)

First, create a membrane with just POPC. With the following command, we're defining a 75x75 A membrane.

```
packmol-memgen --lipids POPC --distxy_fix 75
```

> Note: We are not adding any concentration of ions. Ideally there should be a salt concentration of 0.15M to replicate more accurately the real conditions. However, for the purposes of the tutorial, we're not adding any. Our system is already neutral so it won't necessarily affect the electrostatics.

This process will take some minutes. Once it's done, you can check the resulting membrane with VMD or other visualization software such as [PyMOL](https://www.pymol.org).

```
vmd bilayer_only.pdb
```

#### System topology and transformation to GROMACS

Next, we need to obtain to obtain the Amber force field parameters (version 14SB, ff14SB) for our system. This is done with the processing tool `leap` that will output a .prmtop and .inpcrd files provided of a PDB file.

```
cp ../files/leap.in
```
```
tleap -f leap.in
```

Recent versions of PACKMOL-Memgen tend to produce errors in GROMACS due to structure clashing. Perform a small minimization of the system with the `sander` module to correct it. Input parameters for `sander` are stated in a `system.sander` file. This process can take some time, you can check the progress in the  `system.sanderout` file.

```
../files/system.sander .
```
```
sander -O -i system.sander -o system.sanderout -p system.prmtop -c system.inpcrd -r system.rst -ref system.inpcrd
```

At the moment, we have the required parameters to run the system in [Amber](https://ambermd.org/AmberMD.php) (equally valid). However, we are simulating in [GROMACS](https://manual.gromacs.org/), so we need to transform them to a .top and .gro files. This is done using `amb2gro_top_gro.py`. Additionally, output the minimization end coordinates with the `-b` option.

```
amb2gro_top_gro.py -p system.prmtop -c system.rst -t system_GMX.top -g system_GMX.gro -b system_out.pdb
```
Finally, copy the 'molecular dynamics parameters' files (.mdp) for the following parts.

```
cp -r ../files/mdp .
```

#### Energy minimization

We are going to continue with a short minimization (1000 steps), now in GROMACS. Open the file `mdp/min.mdp` in a text editor to check what you are going to do first.

Execute `gmx grompp` to generate a portable binary run file (.tpr), which contains the starting structure of your simulation, the molecular topology and all the simulation parameters.

```
gmx grompp -f mdp/min.mdp -r system_GMX.gro -c system_GMX.gro -p system_GMX.top -o system_min.tpr -maxwarn 1
```
> Note: Read carefully the NOTES and WARNINGS that are printed during grompp. Usually these are a dead giveaway of wrongdoings that doom your runs to failure.

Then, run `gmx mdrun` to perform the actual calculation.

```
gmx mdrun -deffnm system_min -v
```

#### Equilibration

We are going to perform a three-step equilibration of our system on a NPT ensemble. For this, we need to add a positional restraints statement in the `system_GMX.top` for the lipids which will be gradually tuned down to reach equilibrium. This is done by including the following lines within the lipids' [ moleculetype ] section (named 'system1') but before the water's [ moleculetype ] section begins (named WAT, which can have its own set of positional restraints).

```
#ifdef POSRES
#include "posre.itp"
#endif
```

Thus the `system_GMX.top` file should look something like this.

```
[ moleculetype ]
; Name            nrexcl
system1          3

...

#ifdef POSRES
#include "posre.itp"
#endif

[ moleculetype ]
; Name            nrexcl
WAT          3

...
```

The sequential equilibration runs will be executed through the `equi.sh` shell script. Examine the script, it contains four different sections. Figure out what each of them does. Notice also that each equilibration step has its own `equi*.mdp`.

> Hint: you can use `vimdiff file1 file2` to compare two files and spot the differences more easily.

Once you feel ready, execute the `equi.sh` shell script.

```
chmod +x equi.sh
./equi.sh
```

**IMPORTANT: Without GPU support, the whole equilibration will take a long while. We are not gonna wait for all the steps to finish. Kill the process with `CTRL-C` and copy the `just_popc/*` files from the shared folder provided at the beginning of the class. Then continue with the protocol.**

To visualize the time-evolution of the trajectory of the first equilibration step (equi1) use VMD (or any other visualization software).

```
vmd system_min.gro system_equi1.xtc
```

> Note: The .trr file in GROMACS contains the actual trajectory, but it's best to visualize it using the compressed file .xtc

We are also going to assess some of the variables of this equilibration step. We can do that using GROMACS' analysis tools such as `gmx energy`.

```
gmx energy -f system_equi1.edr -o equi1.xvg
```

Pick the variables we want to assess by typing the following numbers.

`12 14 15 20 21 0`

This way we are selecting the **total energy of the system**, **temperature**, **pressure**, **density** and **volume**. The final zero exits the prompt. To visualize the file you can either use `xmgrace` or execute the python script provided in this repo.

```
python ../files/plot_xvg.py equi1.xvg
```

It outputs a .png image in the same location where the .xvg file is located. See how the different variables change along time until stabilized.

#### Production

We are ready to start producing the proper simulation. We are going to do a 100 ns long unbiased MD simulation. Again, we execute `gmx grompp` and then `gmx mdrun`.

```
gmx grompp -f mdp/prod.mdp -r system_equi3.gro -c system_equi3.gro -p system_GMX.top -o system_prod.tpr -n index.ndx -maxwarn 1
```
> Note: See what has changed between the equi3.mdp and prod.mdp files.
```
gmx mdrun -deffnm system_prod -v
```
**IMPORTANT: Two problems arise here:**
- **1) If you use the provided output for the equilibration `gmx grompp` will likely not go through because the shared system_equi3.gro file and the previously generated system_GMX.top file don't match (inescapable condition). This is true as the two systems were built slightly different.**
- **2) Again, the simulation run takes too long to finish. Proceed with the shared output files.**

#### Analysis

What we're going to do now is a bit of analysis of our membrane. The two typical measurements to assess are the **membrane thickness** and the **area per lipid (APL)**. So as to do that, we're going to use [FATSLiM](http://fatslim.github.io/) a package ready to analyze membrane simulations.

Before doing that you can also check how it visually looks like:

```
vmd system_equi3.gro system_prod.xtc
```

##### Membrane thickness

Before using `fatslim`, we need to define an entry in the index called `headgroups`, indicating the phosphate groups of the lipids. Thus, we need to generate an index with GROMACS first:

```
gmx make_ndx -f system_equi3.gro -o fatslim.ndx<<EOF
a P31
name 8 headgroups
q
EOF
```

This will create an index for fatslim with the atom IDs of the phosphate groups in the POPC molecules. 

To determine membrane thickness the command to be used:

```
fatslim thickness -c system_equi3.gro -t system_prod.xtc -n fatslim.ndx --plot-thickness thickness.xvg
```

The software will give us the thickness per leaflet and for the whole membrane. Moreover, with the option `--plot-thickness` we can obtain a plot of the thickness over time. You can again use `xmgrace` or `plot_xvg.py` to visualize the plot.

##### Area per lipid (APL)

In this case, we will use:

```
fatslim apl -c system_equi3.gro -t system_prod.xtc -n fatslim.ndx --plot-apl apl.xvg
```

And, just like before, the APL per leaflet, for the whole membrane, and a plot over time is generated.

### POPC+CHL

This part of the protocol will be done under `popc+chl/`:

```
cd ../popc+chl
```

Now we're going to simulate a membrane with a ratio of 1 cholesterol molecule per 3 of POPC. Create the membrane with packmol-memgen to begin with:

```
packmol-memgen --lipids POPC:CHL1 --ratio 3:1 --distxy_fix 75
```

Check how it looks like again with `vmd`. Spot where the cholesterol molecules are.

In the same OneDrive link provided before there's also a 100 ns simulation of a similar POPC+CHL membrane system. Repeat the previous analysis to measure the membrane thickness and APL, and compare those to the POPC-only membrane.

```
gmx make_ndx -f system_equi3.gro -o fatslim.ndx<<EOF
r PC | r CHL & a P31 | a O1
name 9 headgroups
q
EOF
```

## Protein-Membrane system

This part of the protocol will be done under `membrane_protein/`:

```
cd ../membrane_protein
```

### Building the system

In here we're going to create a bilayer for a membrane protein. Our membrane protein is a refined structure of the Cannabinoid Receptor 2 (CB2). 

![](files/images/pdbcb2.png)

Run the following command:

```
packmol-memgen --pdb protein.pdb --lipids POPC:CHL1 --ratio 10:1 \
    --dist 12 --dist_wat 15 --salt --salt_c Na+
```

What we're doing now is building a bilayer for our protein, 12 A to the edges of the box in the X and Y axis, 15 A in the Z axis, and with a concentration of salt of 0.15 M (need to specify the cation).

Check the output .PDB with vmd:

```
vmd bilayer_protein.pdb
```

We're not gonna go over again the steps we followed before because it's probably gonna take even longer than the previous parts. So we're going to go straight to the analysis part. 

### Analysis

In the same OneDrive link, a 100 ns trajectory of a system (similar) to the one you've created. That is a CB2 receptor embedded in a 10:1 POPC:CHL membrane, with 0.15 M NaCl. We're going to proceed now to the analysis of some variables concerning the membrane protein.

#### RMSD

We're gonna first measure the RMSD of the C-alpha atoms of our receptor along the trajectory. This can be done with GROMACS:

```
gmx rms -f system_prod.xtc -s system_equi6.gro -o rmsd.xvg
```

Select the "C-alpha" group twice (type "3", press Enter, and type "3" again). GROMACS will automatically align all the coordinates and calculate the RMSD for the C-alpha atoms of our protein.

You can again plot it with the python script.

#### RMSF

Now we're going to calculate the RMSF of the structure througout the simulation. This will help us determine how stable are our transmembrane alpha-helices. 

```
gmx rmsf -f system_prod.xtc -s system_equi6.gro -o rmsf.xvg
```

Select "C-alpha" as well. And again, plot it with whatever you want.

#### Secondary Structure analysis

Finally we're going to perform a simple secondary structure analysis to further assess the stability of the TM helices.

```
gmx dssp -f system_prod.xtc -s system_equi6.gro -o dssp.dat
```

Technically it should know, based on the `index.ndx`, what the atoms of the Protein are. The output `dssp.dat` contains the secondary structure prediction per residue (each character) per frame (each line). Each letter stand for:
- H = α-helix
- B = residue in isolated β-bridge
- E = extended strand, participates in β ladder
- G = 3-helix (310 helix)
- I = 5 helix (π-helix)
- T = hydrogen bonded turn
- S = bend
- ~ = unknown

You can try to plot it yourselves (this won't work with XMGRACE), or use this trashy script:

```
python ../files/plot_dssp.py dssp.dat
```

This will return a `dssp.png` with a heatmap, where each color represents a letter (secondary structure), the sequence in the X-axis and each frame in the Y-axis.

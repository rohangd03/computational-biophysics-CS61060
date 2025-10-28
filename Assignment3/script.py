import sys
sys.path.append("/mnt/d/PyRosetta4.Debug.python311.linux.release-387")
import time
import pyrosetta
pyrosetta.init(extra_options = "-in::file::fullatom -mute all")

from pyrosetta import create_score_function
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.simple_moves import (
    SwitchResidueTypeSetMover,
    ClassicFragmentMover,
)
from pyrosetta.rosetta.core.fragment import ConstantLengthFragSet
from pyrosetta.rosetta.protocols.moves import (
    SequenceMover,
    RepeatMover,
    TrialMover,
    MonteCarlo,
)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta import get_score_function

# ========== STEP 1: CREATE INITIAL POSE ==========
sequence = "MLSDEDFKAFGMTRSAFANLPLWKQQNLKKEKLLF"

ENERGY_FILE = "energy_convergence_data.txt"
N_long_frag_repeats = 1
N_short_frag_repeats = 3
num_refinement_cycles = 5
num_centroid_cycles = 15000
OUTPUT_FILE = "output.pdb"


# ========== PREPARE INITIAL POSE =========================
sequence = "MLSDEDFKAFGMTRSAFANLPLWKQQNLKKEKLLF"
pose = pyrosetta.pose_from_sequence(sequence, "fa_standard")

print(f"Pose created with {pose.total_residue()} residues.")

scorefxn_centroid = create_score_function("score3")
scorefunc_fullatom = get_score_function()

frag9_file = "/home/rohan/aat000_09.frag"
frag3_file = "/home/rohan/aat000_03.frag"

frag9 = ConstantLengthFragSet(9)
frag9.read_fragment_file(frag9_file)
frames9 = frag9.nr_frames()

frag3 = ConstantLengthFragSet(3)
frag3.read_fragment_file(frag3_file)
frames3 = frag3.nr_frames()

if frames9 == 0 or frames3 == 0:
        print("ERROR: Failed to load sufficient fragments from files. Exiting.", file=sys.stderr)
        sys.exit(1)



# ========== STEP 2: LINEARIZE AND CONVERT TO CENTROID ==========
print("Linearizing pose (Setting extended backbone torsional angles) . . . .")

for i in range(1, pose.total_residue() + 1):
    pose.set_phi(i, -150)
    pose.set_psi(i, 150)
    pose.set_omega(i, 180)

print(f"Converting pose to CONTROID MODE . . . ")
switch_to_centroid = SwitchResidueTypeSetMover("centroid")
switch_to_centroid.apply(pose)
print("Pose converted to centroid mode.")



# ========== STEP 3: MOVEMAP AND FRAGMENTS ==========
print(f"Setting up MoveMap . . . .")

movemap = MoveMap()
movemap.set_bb(True)

mover9 = ClassicFragmentMover(frag9, movemap)
mover3 = ClassicFragmentMover(frag3, movemap)


# ========== STEP 4: MONTE CARLO SETUP ==========
ENERGY_FILE = "energy_convergence_data.txt"
N_long_frag_repeats = 1
N_short_frag_repeats = 3
num_refinement_cycles = 5
num_centroid_cycles = 15000
OUTPUT_FILE = "output.pdb"

seq_mover = SequenceMover()
seq_mover.add_mover(RepeatMover(mover9, N_long_frag_repeats))  # 3× 9-mer insertions
seq_mover.add_mover(RepeatMover(mover3, N_long_frag_repeats))  # 5× 3-mer insertions

print(f"Setting up the Monte Carlo simulation parameters . . . ")
kT = 1.0
mc = MonteCarlo(pose, scorefxn_centroid, kT)
trial_mover = TrialMover(seq_mover, mc)

print(f"Running {num_centroid_cycles} cycles of Monte Carlo sampling . . .")
start_ = time.time()
energyValues = []

for i in range(num_centroid_cycles):
      trial_mover.apply(pose)

      lowestScore = mc.lowest_score()
      energyValues.append((i+1, lowestScore))

      if (i + 1) % (num_centroid_cycles // 20 or 1) == 0:
             print(f"Cycle {i+1}/{num_centroid_cycles} . . .  Lowest Score: {lowestScore:.2f} REU")

end_ = time.time()
mc.show_state()

print(f"Centroid sampling completed. The process ran for {end_ - start_} seocnds . . . ")

print(f"Retrieving lowest energy centroid structure . . . ")
mc.recover_low(pose)
final_lowest = mc.lowest_score()
print(f"Lowest centroid energy: {mc.lowest_score():.2f} kcal/mol")


print(f"- - - Saving Energy Convergence Data - - -")
try:
    with open(ENERGY_FILE, "w") as f:
        f.write("Lowest energy data for each cycle \n")
        f.write("Cycle \t Lowest_Energy_REU \n")
        for cycle, energy in energyValues:
            f.write(f"{cycle} \t {energy:.3f} \n")

    print(f"  Energy data saved to '{ENERGY_FILE}'")

except Exception as e:
    print(f"Error: Couldm't write energy data to file")
    print(f"{e}", file=sys.stderr)




# ========== STEP 5: RECOVERY AND CONVERSION ==========
switch_to_fa = SwitchResidueTypeSetMover("fa_standard")
switch_to_fa.apply(pose)
print("Pose converted back to full-atom mode.")

initial_fa_SCORE = scorefunc_fullatom(pose)
print(f"Initial Full-Atom Score (Before Relax): {initial_fa_SCORE:.3f} REU")

# Optional scoring in full-atom (not part of folding)

fa_scorefxn = get_score_function(True)
relax = FastRelax()
relax.set_scorefxn(fa_scorefxn)
relax.apply(pose)


### Refinement - - - Refining using Fast Relax
print(f"Refining using Fast Relax . . .")
relax = FastRelax()
relax.set_scorefxn(scorefunc_fullatom)
relax.max_iter(num_refinement_cycles)


START = time.time()
relax.apply(pose)
END = time.time()

final_fa__score = scorefunc_fullatom(pose)

print(f"FastRelax completed . . . Time for execution: {END - START :.2f} seconds")
print(f"Full-atom energy after conversion: {final_fa__score:.2f}")


# ========== STEP 6: OUTPUT FINAL STRUCTURE ==========

print(f"Savin final structure to '{OUTPUT_FILE}' . . . ")

try:
    pose.dump_pdb(OUTPUT_FILE)
    print(f"Final structure saved!")
except Exception as e:
    print(f"Error saving PDB file: {e}", file = sys.stderr)
    sys.exit(1)



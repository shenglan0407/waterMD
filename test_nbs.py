import matplotlib.pyplot as plt
import mdtraj as md
import water_stats
import numpy as np
import h5py

# from pympler.tracker import SummaryTracker
# tracker = SummaryTracker()

def test_dist(tthd):
    dist=[np.linalg.norm(np.array(tthd[ii])-np.array(tthd[-1])) for ii in range(len(tthd))]
    return dist

traj=md.load_trr('data/nvt-pr_run3.trr',top='data/water-sol_run3.gro')
cut_off=0.5

# reload(water_stats)
w2=water_stats.WaterStats(traj,'run3')
output_file=h5py.File('large_data/run3_nearest_nbs.hdf5','a')
dset_name='cuttoff_%.1f'%cut_off
all_atoms_nbs=[]
for frame in range(1):
    all_atoms_nbs.extend(w2.make_nearest_nbs(cut_off,frame))
print type(all_atoms_nbs[0][0])
output_file.create_dataset('test',data=np.array(all_atoms_nbs))

plt.plot(test_dist(all_atoms_nbs[3]))

plt.show()
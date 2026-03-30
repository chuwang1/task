import matplotlib.image as mpimg
import numpy as np

for file in ["profile_q_p_f_dp_ff_comparison.png", "comparison_overview_big.png"]:
    try:
        img = mpimg.imread(file)
        print(f"{file} shape: {img.shape}")
        # check if it's completely white/blank
        if np.all(img[:,:,:3] == 1.0):
            print(f"  -> {file} is completely white!")
        else:
            std_dev = np.std(img[:,:,:3])
            print(f"  -> {file} has content. Std dev: {std_dev}")
    except Exception as e:
        print(f"Error loading {file}: {e}")

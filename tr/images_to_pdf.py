import os
import glob
from PIL import Image

def create_pdf(output_filename):
    # Find all PNGs starting with tr_data_
    images = sorted(glob.glob("tr_data_*.png"))
    
    if not images:
        print("No images found.")
        return

    print(f"Found {len(images)} images. Combining into {output_filename}...")
    
    pil_images = []
    for img_path in images:
        try:
            img = Image.open(img_path)
            # Convert to RGB (PDF doesn't support RGBA with transparency well in this simplified mode)
            if img.mode == 'RGBA':
                img = img.convert('RGB')
            pil_images.append(img)
        except Exception as e:
            print(f"Failed to load {img_path}: {e}")

    if pil_images:
        # Save the first image and append the rest
        first_image = pil_images[0]
        if len(pil_images) > 1:
            first_image.save(output_filename, "PDF", resolution=100.0, save_all=True, append_images=pil_images[1:])
        else:
            first_image.save(output_filename, "PDF", resolution=100.0)
        print(f"Successfully created {output_filename}")
    else:
        print("No valid images to save.")

if __name__ == "__main__":
    create_pdf("tr_simulation_results.pdf")

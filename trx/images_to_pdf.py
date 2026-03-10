import os
import glob
import math
import sys
from PIL import Image

def get_grid_dims(n):
    """
    Calculate generic grid dimensions (rows, cols) for n items
    trying to keep it somewhat rectangular or square.
    """
    if n <= 1: return 1, 1
    if n == 2: return 1, 2 # 1 row, 2 cols (side by side)
    if n == 3: return 1, 3 # 1 row, 3 cols
    if n == 4: return 2, 2
    if n <= 6: return 2, 3 # 2 rows, 3 cols
    if n <= 8: return 2, 4
    if n <= 9: return 3, 3
    if n <= 12: return 3, 4
    
    # Fallback generic: square-ish
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    return rows, cols

def combine_images_on_page(images, rows, cols, bg_color=(255, 255, 255)):
    if not images: return None
    
    # Assume all images are roughly same size or take max
    # If images are vastly different, this generic grid might look weird, 
    # but for loop plots they are usually identical.
    widths = [img.width for img in images]
    heights = [img.height for img in images]
    max_w = max(widths)
    max_h = max(heights)
    
    # Total canvas size
    page_w = max_w * cols
    page_h = max_h * rows
    
    page_img = Image.new('RGB', (page_w, page_h), bg_color)
    
    for idx, img in enumerate(images):
        r = idx // cols
        c = idx % cols
        # Calculate position
        x = c * max_w
        y = r * max_h
        
        # Center image in cell
        offset_x = x + (max_w - img.width) // 2
        offset_y = y + (max_h - img.height) // 2
        
        page_img.paste(img, (offset_x, offset_y))
        
    return page_img

def create_pdf(output_filename):
    # Find all PNGs starting with tr_data_
    # Sort by number in filename if possible (tr_data_1.png vs tr_data_10.png)
    # Default string sort tr_data_1, tr_data_10, tr_data_2... is usually not what we want.
    
    image_files = glob.glob("tr_data_*.png")
    if not image_files:
        print("No images found.")
        return

    # Try to extract number for sorting
    def sort_key(chk):
        # Extract numbers: tr_data_123.png -> 123
        base = os.path.basename(chk)
        name, ext = os.path.splitext(base)
        parts = name.split('_')
        for p in parts:
            if p.isdigit():
                return int(p)
        return chk
    
    images = sorted(image_files, key=sort_key)

    print(f"Found {len(images)} images.")
    
    # User input
    # Check if arg is passed via command line, else ask input
    if len(sys.argv) > 1 and sys.argv[1].isdigit():
         num_per_page = int(sys.argv[1])
         print(f"Using command line argument: {num_per_page} images per page.")
    else:
        try:
            print("Entrer number of images per page (default 6):")
            # Flush stdout to ensure prompt appears before input
            sys.stdout.flush()
            n_str = sys.stdin.readline().strip() 
            # Note: in some envs input() is tricky, sys.stdin is safer if piped
            if not n_str:
                num_per_page = 6
            else:
                num_per_page = int(n_str)
        except ValueError:
            print("Invalid input, defaulting to 6.")
            num_per_page = 6

    print(f"Combining {len(images)} images into {output_filename} with {num_per_page} per page...")

    chunk_size = num_per_page
    rows, cols = get_grid_dims(num_per_page)
    print(f"Grid Layout: {rows} rows x {cols} cols")

    pil_pages = []
    
    for i in range(0, len(images), chunk_size):
        chunk_paths = images[i:i + chunk_size]
        chunk_imgs = []
        for p in chunk_paths:
            try:
                img = Image.open(p)
                if img.mode == 'RGBA':
                    img = img.convert('RGB')
                chunk_imgs.append(img)
            except Exception as e:
                print(f"Error reading {p}: {e}")
        
        if chunk_imgs:
            page = combine_images_on_page(chunk_imgs, rows, cols)
            pil_pages.append(page)
    
    if pil_pages:
        pil_pages[0].save(output_filename, "PDF", resolution=100.0, save_all=True, append_images=pil_pages[1:])
        print(f"Successfully created {output_filename}")
    else:
        print("No pages created.")

if __name__ == "__main__":
    create_pdf("tr_simulation_results.pdf")

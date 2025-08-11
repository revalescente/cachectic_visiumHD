import os
import re
import json
import numpy as np
from tqdm import tqdm
from skimage import io

def load_qupath_settings(json_path):
    """Carica le impostazioni dei canali dal file JSON di QuPath."""
    with open(json_path, 'r') as f:
        settings = json.load(f)
    # Filtra solo i canali attivi (isShowing=True) e ordinati come R, G, B
    rgb_channels = [
        ch for ch in settings['channels'] 
        if ch['isShowing'] and ch['name'] in ['Red', 'Green', 'Blue']
    ]
    rgb_channels.sort(key=lambda x: ['Red', 'Green', 'Blue'].index(x['name']))
    return rgb_channels

def apply_qupath_rgb(image, channels):
    """Applica le impostazioni di colore di QuPath a un'immagine."""
    normalized = np.zeros_like(image[..., :3], dtype=np.float32)  # Ignora il 4° canale
    for i, ch in enumerate(channels):
        min_val, max_val = ch['minDisplay'], ch['maxDisplay']
        channel_data = image[..., i].astype(np.float32)
        normalized[..., i] = np.clip((channel_data - min_val) / (max_val - min_val), 0, 1)
    
    rgb_output = np.zeros_like(normalized)
    for i, ch in enumerate(channels):
        color = np.array([ch['color']['red'], ch['color']['green'], ch['color']['blue']]) / 255.0
        rgb_output += normalized[..., i:i+1] * color
    return np.clip(rgb_output, 0, 1)

def process_images_with_custom_jsons(image_json_mapping, output_dir):
    """
    Elabora immagini applicando a ciascuna un JSON specifico.
    Args:
        image_json_mapping (dict): Dizionario {percorso_immagine: percorso_json}
        output_dir (str): Cartella di output
    """
    os.makedirs(output_dir, exist_ok=True)
    

    for image_path, json_path in tqdm(image_json_mapping.items(), desc="Work in progress"):
        # Carica immagine e impostazioni JSON
        image = io.imread(image_path)
        channels = load_qupath_settings(json_path)
        
        # Applica regole di colore
        rgb_image = apply_qupath_rgb(image, channels)
        
        # Estrai numero blocco dal nome file
        match = re.search(r'Blocco(\d)_20x\.tif', os.path.basename(image_path))
        if match:
            blocco_num = match.group(1)
            output_filename = f"pp_blocco{blocco_num}_20x.tif"
            output_path = os.path.join(output_dir, output_filename)
            
            # Salva immagine
            io.imsave(output_path, (rgb_image * 255).astype(np.uint8))
        else:
            print(f"Pattern non trovato per file: {image_path}")

# --- ESEMPIO DI USO ---
if __name__ == "__main__":
    # Mappatura: {immagine: json_da_applicare}
    image_json_mapping = {
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco1_20x.tif": "data/visual_json/cachetic_correct_visual.json",
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco2_20x.tif": "data/visual_json/cachetic_correct_visual.json",
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco3_20x.tif": "data/visual_json/cachetic_correct_visual.json",
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco4_20x.tif": "data/visual_json/cachetic_correct_visual.json",
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco5_20x.tif": "data/visual_json/cachetic_correct_visual.json",
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco6_20x.tif": "data/visual_json/cachetic_correct_visual.json",
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco7_20x.tif": "data/visual_json/cachetic_correct_visual_b7.json",
        "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco9_20x.tif": "data/visual_json/cachetic_correct_visual_b7.json",
    }
    output_dir = "HE_images/color_corrected"
    
    process_images_with_custom_jsons(image_json_mapping, output_dir)
    print("Work done!")

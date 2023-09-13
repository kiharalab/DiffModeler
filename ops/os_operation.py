import os
def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        print (path+" created")
        os.makedirs(path)
        return True
    else:
        print (path+' existed')
        return False
def execCmd(cmd):
    r = os.popen(cmd)
    text = r.read()
    r.close()
    return text

def run_command(cmd):
    os.system(cmd)

import shutil

def copy_directory(source_dir, destination_dir):
    """
    Copies a directory and its contents to the specified destination directory.

    Args:
        source_dir (str): Path to the source directory.
        destination_dir (str): Path to the destination directory.

    Returns:
        None
    """
    try:
        shutil.copytree(source_dir, destination_dir)
        print(f"Directory {source_dir} copied to {destination_dir}")
    except shutil.Error as e:
        print(f"Error copying directory: {e}")
    except OSError as e:
        print(f"Error creating destination directory: {e}")
    return destination_dir
def functToDeleteItems(fullPathToDir):
   for itemsInDir in os.listdir(fullPathToDir):
        if os.path.isdir(os.path.join(fullPathToDir, itemsInDir)):
            functToDeleteItems(os.path.join(fullPathToDir, itemsInDir))
        else:
            os.remove(os.path.join(fullPathToDir,itemsInDir))
import gzip
def unzip_gz(file_path):
    new_path = file_path.replace(".gz","")
    g_file = gzip.GzipFile(file_path)
    with open(new_path,"wb+") as file:
        file.write(g_file.read())
    g_file.close()
    return new_path
import os
import tarfile
import zipfile
import gzip

def extract_compressed_file(file_path, extract_dir):
    """
    Extracts a compressed file to the specified directory.

    Args:
        file_path (str): Path to the compressed file.
        extract_dir (str): Directory where the contents will be extracted.

    Returns:
        None
    """
    base_name = os.path.basename(file_path)

    if file_path.endswith('.tar.gz') or file_path.endswith('.tgz'):
        with tarfile.open(file_path) as tar:
            tar.extractall(path=extract_dir)
    elif file_path.endswith('.tar'):
        with tarfile.open(file_path) as tar:
            tar.extractall(path=extract_dir)
    elif file_path.endswith('.gz'):
        output_path = os.path.join(extract_dir, os.path.splitext(base_name)[0])
        with open(output_path, 'wb') as out_f, gzip.open(file_path, 'rb') as gz_f:
            out_f.write(gz_f.read())
    elif file_path.endswith('.zip'):
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
    else:
        print("Unsupported zipped file format %s"%file_path)
        return

    print(f"Successfully extracted {file_path} to {extract_dir}")

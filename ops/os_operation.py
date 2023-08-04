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

import gzip
def unzip_gz(file_path):
    new_path = file_path.replace(".gz","")
    g_file = gzip.GzipFile(file_path)
    with open(new_path,"wb+") as file:
        file.write(g_file.read())
    g_file.close()
    return new_path

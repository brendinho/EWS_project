# copied from stackexchange
# https://stackoverflow.com/questions/748675/finding-duplicate-files-and-removing-them
# user Todor
# last edited on SE: 17 Oct 2017
# taken to Python3 and edited by Brendon Phillips: 02 Jan 2018

#!/usr/bin/env python3
import sys
import os
import hashlib


def chunk_reader(fobj, chunk_size=1024):
    """Generator that reads a file in chunks of bytes"""
    while True:
        chunk = fobj.read(chunk_size)
        if not chunk:
            return
        yield chunk


def get_hash(filename, first_chunk_only=False, hash=hashlib.sha1):
    hashobj = hash()
    file_object = open(filename, 'rb')

    if first_chunk_only:
        hashobj.update(file_object.read(1024))
    else:
        for chunk in chunk_reader(file_object):
            hashobj.update(chunk)
    hashed = hashobj.digest()

    file_object.close()
    return hashed


# def check_for_duplicates(paths, hash=hashlib.sha1):
def check_for_duplicates(path, hash=hashlib.sha1):
    hashes_by_size = {}
    hashes_on_1k = {}
    hashes_full = {}

    # folder_for_duplicate_files = "{0}/Duplicate_Files".format(path);
	#
    # if not os.path.exists(folder_for_duplicate_files):
    #     os.makedirs(folder_for_duplicate_files)

    for filename in os.listdir(path):

        full_path = "{0}/{1}".format(path, filename);

        # if not filename.endswith(".bin"):
        #     continue
        try:
            file_size = os.path.getsize(full_path);
        except (OSError,):
            continue

        if hashes_by_size.get(file_size):
            hashes_by_size[file_size].append(full_path)
        else:
            hashes_by_size[file_size] = []  # create the list for this file size
            hashes_by_size[file_size].append(full_path)

    # For all files with the same file size, get their hash on the 1st 1024 bytes
    for __, files in hashes_by_size.items():
        if len(files) < 2:
            continue    # this file size is unique, no need to spend cpy cycles on it

        for filename in files:
            small_hash = get_hash(filename, first_chunk_only=True)

            if hashes_on_1k.get(small_hash):
                hashes_on_1k[small_hash].append(filename)
            else:
                hashes_on_1k[small_hash] = []          # create the list for this 1k hash
                hashes_on_1k[small_hash].append(filename)

    # For all files with the hash on the 1st 1024 bytes, get their hash on the full file - collisions will be duplicates
    for __, files in hashes_on_1k.items():
        if len(files) < 2:
            continue    # this hash of first 1k file bytes is unique, no need to spend cpy cycles on it

        for filename in files:
            full_hash = get_hash(filename, first_chunk_only=False)
            already_existing_file = hashes_full.get(full_hash)

            if already_existing_file:
                # print( "Duplicate found: {0} and {1}".format(filename, already_existing_file) );
                # os.rename( filename, "{0}/{1}".format(folder_for_duplicate_files, filename.split("/")[-1] ));
                os.remove(filename)
            else:
                hashes_full[full_hash] = filename

# if sys.argv[1:]:
#     check_for_duplicates(sys.argv[1:])
# else:
#     print "Please pass the paths to check as parameters to the script"

# check_for_duplicates("/mnt/c/Users/brendon/Desktop/All_Data");
# check_for_duplicates("/Users/b2philli/Desktop/Frames")
# check_for_duplicates("/home/b2philli/Desktop/Files/");
# check_for_duplicates(sys.argv[1])
check_for_duplicates("/media/b2philli/Simulations/Hesitance")

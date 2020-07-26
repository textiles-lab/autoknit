#!/usr/bin/env python3

import sys

import urllib.request
import os
import subprocess
import shutil
import platform
import re

print("Fetching jam...")

work_folder = "work"

if not os.path.exists(work_folder):
	print("Creating work folder '" + work_folder + "'")
	os.mkdir(work_folder)

def run_command(args,cwd=None,env=None):
	print("  Running `\"" + '" "'.join(args) + "\"`")
	subprocess.run(args,cwd=cwd,check=True,env=env)

def unzip_file(filename, folder):
	run_command([
		"C:\\Program Files\\7-Zip\\7z.exe",
		"x",
		"-o" + folder,
		filename
	])

def fetch_file(url, filename, checksum=None):
	if os.path.exists(filename):
		print("  File '" + filename + "' exists; TODO: check checksum.")
		return
	print("  Fetching '" + url + "' => '" + filename + "'")
	urllib.request.urlretrieve(url, filename)


jam_file = 'ftjam-2.5.2-win32.zip'
jam_url = 'https://sourceforge.net/projects/freetype/files/ftjam/2.5.2/' + jam_file + '/download'

def fetch_jam():
	fetch_file(jam_url, work_folder + "/" + jam_file)
	unzip_file(work_folder + "/" + jam_file, work_folder)

fetch_jam()

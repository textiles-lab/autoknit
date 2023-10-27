#!/usr/bin/env python3

import sys
import urllib.request
urllib.request.urlretrieve(sys.argv[1], sys.argv[2])

import os

def assert_file_exists(fn,test_nonempty=False):
	assert os.path.isfile(fn), "File '{}' does not exist.".format(fn)

def assert_file_nonempty(fn):
	statinfo = os.stat(fn)
	file_size = statinfo.st_size
	assert file_size > 0, "File '{}' is empty".format(fn)

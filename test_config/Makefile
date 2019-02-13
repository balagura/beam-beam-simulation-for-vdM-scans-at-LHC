.PHONY : test
test : test_config
	test_config test_config.txt

test_config : ../config.hh ../config.C test_config.C
	g++ ../config.C test_config.C -I.. -o test_config

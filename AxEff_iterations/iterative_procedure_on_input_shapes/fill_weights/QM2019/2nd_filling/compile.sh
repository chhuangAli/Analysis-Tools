SRC=outputAna.cpp
g++ -o GOrun $SRC `root-config --cflags --libs` -I/home/chunlu/local/alice_2017/sw/ubuntu1604_x86-64/AliRoot/0_ROOT6-1/include -L/home/chunlu/local/alice_2017/sw/ubuntu1604_x86-64/AliRoot/0_ROOT6-1/lib -lSTEERBase

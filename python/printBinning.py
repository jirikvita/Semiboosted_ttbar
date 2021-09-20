#!/usr/bin/python

a = [ 300 + 40.*i for i in range(0, 40)]
s = ''
for b in a:
    s = s + ', {}'.format(b)
print s

a = [ 25.*i for i in range(1, 50)]
s = ''
for b in a:
    s = s + ', {}'.format(b)
print s



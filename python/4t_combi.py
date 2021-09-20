#!/usr/bin/python

# jk 3.12.2018


from __future__ import print_function

# resolved, semiboosted, boosted
regimes = ['r', 's', 'b']

counts = {}

for i in regimes:
    for j in regimes:
        for k in regimes:
            for l in regimes:
                regime = i+j+k+l
                #print('"{:}"'.format(regime,))
                regime = ''.join(sorted(regime)) # trick needed!
                print('Working on "{:}"'.format(regime,))
                            
                try:
                    counts[regime] = counts[regime] + 1
                except:
                    counts[regime] = 1

print('OK, some 4t final state combinatorics based on r=resolved, s=semiboosted and b=boosted combinatorics ;-)')
print('Resulting unique combinations')
#print(counts)
sum = 0
for count in counts:
    val = counts[count]
    print('{:} : {:}'.format(count, val))
    sum = sum + val
print('Unique combinations: {}'.format(len(counts)),)
print('All combinations: {}'.format(sum,))


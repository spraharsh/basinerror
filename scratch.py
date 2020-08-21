import numpy as np

point = [0, 0]
blah = 0.1
nmesh = 1
a = np.linspace(point[0] - blah/2, point[0]+ blah/2, 1)

blah = np.arange(0, 18, 1)

b = blah.reshape(18//3, 3)


bn = map(lambda x: x-np.array([0, 1, 2]), b)


print(blah)
blah = np.array(list(bn))
print(blah.shape)
print(blah)

print(blah%4)
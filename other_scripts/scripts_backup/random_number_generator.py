import random

a = list(range(1, 1001))

abundance_1 = random.sample(a, 30)
abundance_2 = random.sample(a, 30)
abundance_3 = random.sample(a, 30)

print('%s\t%s' % (sum(abundance_1), abundance_1))
print('%s\t%s' % (sum(abundance_2), abundance_2))
print('%s\t%s' % (sum(abundance_3), abundance_3))


genomes = ['AAM.fna', 'AKV.fna', 'AMAC.fna', 'AMAU.fna', 'AMS.fna', 'ARL.fna', 'ARS.fna', 'ASJ.fna', 'ASN.fna', 'ATM.fna', 'BAD.fna', 'BDS.fna', 'BGC.fna', 'BHS.fna', 'BNM.fna', 'BRT.fna', 'BSA.fna', 'BSD.fna', 'BSL.fna', 'BTK.fna', 'GAV.fna', 'GCR.fna', 'GGA.fna', 'GGN.fna', 'GHN.fna', 'GPS.fna', 'GRO.fna', 'GSM.fna', 'GSO.fna', 'GXN.fna']
abundance_a = [145, 489, 464, 461, 284, 163, 498, 216, 188, 273, 346, 295, 285, 361, 261, 67, 364, 219, 482, 234, 348, 282, 441, 12, 442, 77, 42, 481, 254, 197]
abundance_b = [410, 230, 5, 265, 170, 492, 380, 120, 247, 375, 248, 31, 108, 188, 1, 251, 305, 77, 363, 373, 229, 415, 29, 425, 364, 487, 245, 25, 110, 209]
abundance_c = [86, 12, 166, 183, 58, 233, 434, 206, 328, 159, 145, 127, 232, 5, 2, 373, 93, 364, 124, 72, 200, 162, 106, 224, 367, 143, 286, 175, 304, 396]

n = 0
for each in genomes:
    print('%s\t%s' % (each, abundance_c[n]))
    n += 1


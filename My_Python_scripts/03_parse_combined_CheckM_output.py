
for each in open('/Users/songweizhi/Desktop/64genomes_95_completeness_quality.txt'):

    if not ((each.startswith('---')) or (each.startswith('  Bin Id'))):
        print(each.strip())

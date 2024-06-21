import os
content = {}
with open('/dssg/home/acct-trench/trench/DATA/MEER_RAW/metadata.tsv', 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        content[line[0]] = content.get(line[0], []) + [line[1]]

samples = {}
with open('tst_low_depth/data/list', 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        samples[line[0]] = content[line[0]]

# write to .path fike
for key, value in samples.items():
    with open('tst_low_depth/data/samples/' + key + '.path', 'w') as f:
        for line in value:
            f.write('{0}\n'.format(line))
    print('{0} created'.format(key))

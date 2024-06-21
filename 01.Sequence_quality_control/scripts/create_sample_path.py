import os
content = {}
with open('metadata.tsv', 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        content[line[0]] = content.get(line[0], []) + [line[1]]

samples = {}
# list is a text file contains a list of sample names to be included
with open('data/list', 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        samples[line[0]] = content[line[0]]

# write to .path fike
for key, value in samples.items():
    with open('data/samples/' + key + '.path', 'w') as f:
        for line in value:
            f.write('{0}\n'.format(line))
    print('{0} created'.format(key))

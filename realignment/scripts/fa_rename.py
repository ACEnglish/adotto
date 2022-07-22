import sys

new_name = sys.argv[1]
old_name = sys.stdin.readline()
print(f'>{new_name}')
for line in sys.stdin:
    sys.stdout.write(line)

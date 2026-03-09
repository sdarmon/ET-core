import sys

Arg = sys.argv[:]

def homomorphic_compression(read,w):
    #Removing consecutive duplicates window by window
    S = read[0:w]
    current_w = read[0:w]
    for i in range(1, len(read)-w+1):
        next_w = read[i:i+w]
        if current_w != next_w:
            S += read[i+w-1]
            current_w = next_w
    return S

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " reads.fq reads_out.fq window")
    exit()

cycle = 0
w = int(Arg[3])
if len(Arg) == 4:
    with open(Arg[1], 'r') as f:
        with open(Arg[2], 'w') as f_out:    
            for line in f:
                cycle = (cycle + 1) % 4
                if len(line) != "" :
                    if cycle != 2:
                        f_out.write(line)
                        continue
                    else:
                        f_out.write(homomorphic_compression(line,w))
#!/usr/bin/awk -f

/^>/ {
    # header line
    id = substr($1, 2)  # remove '>'
    for (i = 2; i <= NF; i++) {
        # look for fields starting with L:
        if ($i ~ /^L:/) {
            n = split($i, a, ":")  # a[1]=L, a[2]=sign1, a[3]=number, a[4]=sign2
            s1 = a[2]
            num = a[3]
            s2 = a[4]

            if (s1 == "+" && s2 == "+")      type = "FF"
            else if (s1 == "+" && s2 == "-") type = "FR"
            else if (s1 == "-" && s2 == "+") type = "RF"
            else if (s1 == "-" && s2 == "-") type = "RR"
            else                              type = "??"

            printf("%s\t%s\t%s\n", id, num, type)
        }
    }
}


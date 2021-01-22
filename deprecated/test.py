
bins=[]
myfile="Z:/jobs/Qfep_NEW/FEP1_md_1000_0000/FEP1_md_1000_0000.en"
with open(myfile, "rb") as f:
    byte = f.read(1)
    while byte !=b"":
        # Do stuff with byte.
        byte = f.read(1)
        bins.append(byte)
print(bins)
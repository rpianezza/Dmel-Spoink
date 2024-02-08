import os
import re
import textwrap
import subprocess
import pandas as pd


tsv_mel = '/Users/ascarpa/Desktop/megatest/dmel.tsv'
tsv_wil = '/Users/ascarpa/Desktop/megatest/dwil.tsv'

df_mel = pd.read_csv(tsv_mel, sep='\t')
df_wil = pd.read_csv(tsv_wil, sep='\t')

folder_mel = '/Users/ascarpa/Desktop/megatest/Spoink/D.mel.Es_Ten/run_diptera_odb10/busco_sequences/single_copy_busco_sequences/'
folder_wil = '/Users/ascarpa/Desktop/megatest/Spoink/D.willistoni.17/run_diptera_odb10/busco_sequences/single_copy_busco_sequences/'

dSdN_file = '/Users/ascarpa/Desktop/megatest/seqfaa/dn_ds.txt'

n_ort = 0

for value_mel in df_mel.iloc[:, 0]:
    if value_mel in df_wil.iloc[:, 0].values:
        row_position = df_wil.index[df_wil.iloc[:, 0] == value_mel].tolist()[0]
        
        seq_name = df_wil.iloc[row_position, 0]
        seq_faa = f"{seq_name}.faa"
        seq_fna = f"{seq_name}.fna"
        faa_mel_wil = f"/Users/ascarpa/Desktop/megatest/seqfaa/{seq_name}_2.faa"
        fna_mel_wil = f"/Users/ascarpa/Desktop/megatest/seqfaa/{seq_name}_2.fna"

        # Concatenate the protein sequences
        if os.path.exists(os.path.join(folder_mel, seq_faa)) and os.path.exists(os.path.join(folder_wil, seq_faa)):
            with open(os.path.join(folder_mel, seq_faa), 'r') as file1, open(os.path.join(folder_wil, seq_faa), 'r') as file2:
                seq_aa_mel = file1.read()
                seq_aa_wil = file2.read()
            
            with open(faa_mel_wil, 'w') as output_f:
                output_f.write(seq_aa_mel)
                output_f.write(seq_aa_wil)

        else:
            print(f"Files not found for sequence: {seq_name}")


        # Concatenate the DNA sequences
        if os.path.exists(os.path.join(folder_mel, seq_fna)) and os.path.exists(os.path.join(folder_wil, seq_fna)):
            with open(os.path.join(folder_mel, seq_fna), 'r') as file1, open(os.path.join(folder_wil, seq_fna), 'r') as file2:
                seq_DNA_mel = file1.read()
                seq_DNA_wil = file2.read()
            
            with open(fna_mel_wil, 'w') as output_f:
                output_f.write(seq_DNA_mel)
                output_f.write(seq_DNA_wil)
            n_ort += 1
        else:
            print(f"Files not found for sequence: {seq_name}")

        # Align the protein sequences   
        faa_mel_wil_clu = faa_mel_wil + "_clu"
        clustalo_command = ["clustalo", "-i", faa_mel_wil, "-o", faa_mel_wil_clu]
        try:
            subprocess.run(clustalo_command, check=True)
            print("clustalo completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running clustalo: {e}")

        # Align the codons
        fna = os.path.join("/Users/ascarpa/Desktop/megatest/seqfaa/", seq_faa)
        fna_CDS = os.path.splitext(fna)[0] + "_CDS.fna"
        pal2nal_command = ["/Users/ascarpa/Desktop/megatest/pal2nal.v14/pal2nal.pl", faa_mel_wil_clu, fna_mel_wil, "-output", "fasta"]

        try:
            pal2nal_output = subprocess.run(pal2nal_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
            
            with open(fna_CDS, 'w') as output_file:
                output_file.write(pal2nal_output)
            print("pal2nal completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running pal2nal: {e}")


        # Write ctl file
        file_content = textwrap.dedent(f"""
        seqfile = {fna_CDS}
        outfile = {seq_name}
        noisy = 0
        verbose = 1
        runmode = -2
        seqtype = 1
        CodonFreq = 2
        ndata = 10
        clock = 0
        aaDist = 0
        aaRatefile = /Users/ascarpa/Desktop/megatest/paml-old/dat/jones.dat
        model = 2
        NSsites = 0
        icode = 0
        Mgene = 0
        fix_kappa = 0
        kappa = 2
        fix_omega = 0
        omega = .4
        fix_alpha = 1
        alpha = 0
        Malpha = 0
        ncatG = 8
        getSE = 0
        RateAncestor = 1
        Small_Diff = .5e-6
        cleandata = 1
        *fix_blength = -1
        method = 0
        """)

        ctl = f"/Users/ascarpa/Desktop/megatest/seqfaa/ctl/{seq_name}.ctl"
        with open(ctl, 'w') as file:
            file.write(file_content)


        # paml
        output_paml = '/Users/ascarpa/Desktop/megatest/seqfaa/paml'
        paml_command = ["/Users/ascarpa/Desktop/megatest/paml-old/codeml", ctl]

        try:
            subprocess.run(paml_command, check=True, cwd=output_paml)
            print("Command executed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running the command: {e}")


        # Extract the values for dN/dS, dN and dS
        seq_final = '/Users/ascarpa/Desktop/megatest/seqfaa/paml/' + seq_name
        with open(seq_final, 'r') as file:
            dNdS_line = file.readlines()[-4]

        dN_dS_pattern = re.compile(r'dN/dS=  (\d+\.\d+)')
        dN_pattern = re.compile(r'dN = (\d+\.\d+)')
        dS_pattern = re.compile(r'dS = (\d+\.\d+)')

        dN_dS_match = dN_dS_pattern.search(dNdS_line)
        dN_match = dN_pattern.search(dNdS_line)
        dS_match = dS_pattern.search(dNdS_line)

        if dN_dS_match:
            dN_dS_value = float(dN_dS_match.group(1))
        if dN_match:
            dN_value = float(dN_match.group(1))
        if dS_match:
            dS_value = float(dS_match.group(1))


        if os.path.exists(dSdN_file):
            with open(dSdN_file, 'a') as file:
                file.write(f"{dN_dS_value}\t{dN_value}\t{dS_value}\n")
        else:
            with open(dSdN_file, 'w') as file:
                file.write("dN/dS\t dN\t dS\n")
                file.write(f"{dN_dS_value}\t{dN_value}\t{dS_value}\n")

print('Number of orthologous genes', n_ort)
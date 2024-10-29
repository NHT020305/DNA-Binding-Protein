# Dataset Description
There are some proteins that interact with DNA and controls many characteristics of the DNA. The impact is so powerful that it even affects drug design at industry level. Our goal in this project is to devise a computational method that is able to distinguish DNA binding proteins (positive class) from normal proteins (negative class). There are protein structure based laboratory methods for this, but they are time consuming, expensive and worst of all not even accurate. What we want to do is: design a machine learning based approach to do this task given the protein sequences in a cheap and easy way. You do not need any biology knowledge to work on this problem. But if you want to explore protein biology and relevant feature extraction tools, you are most welcome to do so.

# Dataset Details
### Two files - (1) Train.fasta and (2) Test.fasta
Each sequence consists of a combination of 20 unique characters (each character represents one unique amino acid)
##### Example: 
"MLTPRDENNEGDAMPMLKKPRYSSLSGQSTNITYQEHTISREERAAAVGRHEGFRGCTIWFTGLSGAGKTTISFALERTLNKLGIPCYGLDGDNIRHGLCKNLGFSKEDRQENIRRVAEVAKLFADSGMICLAAFISPFQEDRLDARKIHESENVKFIEVHVSTTLEVCEQRDPKPSELYKKARAGQILGFTGIDSAYEPPENAEIILDAGKDGVQQCVQKVLDHLESKGLLPEQIPDVPAVRELFVSDDLTVAELLKESQNLPTVELTKVDLQWLQVLAEGWATPLSGFMRERQYLQSMHFGQLLDLKHKVAFVGEKSDDKEDSWPMMDDINQSIPIVLPISDDVKKGLEGVTRIALKYNGQVYAILSDPEIFEHRKDERVCRQFGTNDPRHPAVAQVLESGNWLLGGDVAVVQKIQFNDGLDKYRKTPNELRAIFAEKNADAVFAFQLRNPIHNGHALLMRDTREKLLAEHKNPILLLHPLGGWTKDDDVPLDIRIKQHEAVIAERVLDPEWTVLSIFPSPMMYAGPTEVQWHARSRIAAGIQHYIVGRDPAGIQKPGSPDALYETTHGAKVLSMAPGLSALHILPFRVAAYDKTAKKMSFFDTSRKEDFENISGTKMRGLARNGDTPPEGFMAPTAWEVLAGYYKSLQNSN"
### The class label is provided in the header of each protein sequence in the fasta files as "label_1" and "label_0" for positive and negative sequences, respectively
##### Example: ">seq_17848_label_1" (positive sequence) ">seq_17860_label_0" (negative sequence)
### A lot more negative sequences compared to positive sequences in both train and test set
### The protein sequences are of variable length (50 to 3000 characters)

# Expected Task Description
### You need to train and tune your model using Train.fasta
### Finally, you need to test on the Test.fasta
### As test performance metric, you need use sensitivity, specificity and MCC score
### Remember to explore:
##### character level language models
##### various feature extraction techniques
##### class balancing methods during training

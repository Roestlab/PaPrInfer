import sys
from pyopenms import IdXMLFile


def main(epifany_file):
    prot_ids = []
    pep_ids = []

    IdXMLFile().load(epifany_file, prot_ids,
                     pep_ids)



if __name__ == "__main__":
    main(sys.argv[1])

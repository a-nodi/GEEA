"""
__auther__ = https://github.com/a-nodi
Icon made by Freepik from www.flaticon.com (resource/icon.png)

nsa is Nucleic Sequence Analyze
snct is Specific Nucleo base Combine and Translation
rnst is Repeated Nucleic Sequence Translation

"""
import os
import sys
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QSize, pyqtSlot
from loader import load_setting
from analyser import *

# load setting and resource
dict_of_setting = load_setting()
dict_of_resource_dir = load_resource_dir()

CODON_TABLE_DIR = dict_of_resource_dir['codon_table.json']
CODON_TABLE = parse_json(CODON_TABLE_DIR)

SNCT_CODON_TABLE_DIR = dict_of_resource_dir['snct_codon_table.json']
SNCT_CODON_TABLE = parse_json(SNCT_CODON_TABLE_DIR)

AMINO_ACID_ACRONYMS_TO_FULLNAME_TABLE_DIR = dict_of_resource_dir['amino_acid_acronyms_to_fullname.json']
AMINO_ACID_ACRONYMS_TO_FULLNAME_TABLE = parse_json(AMINO_ACID_ACRONYMS_TO_FULLNAME_TABLE_DIR)

ICON_DIR = dict_of_resource_dir['icon_dir']

ui_dir = os.path.join(os.getcwd(), dict_of_resource_dir['ui_dir'])
form_class = uic.loadUiType(ui_dir)[0]  # load Qt UI file

# defining constant that used in gene analysing
LINE_EDIT_1 = 1
LINE_EDIT_2 = 2
FIVEEND_TO_THIRDEND = '53'
THIRDEND_TO_FIVEEND = '35'
IS_DNA = True

FIRST_BASE = 1
SECOND_BASE = 2
THIRD_BASE = 3

U_INDEX = 0
C_INDEX = 1
A_INDEX = 2
G_INDEX = 3


class MainWindow(QMainWindow, form_class):
    """
    GUI class
    """
    def __init__(self):

        super().__init__()

        self.nsa_input1_is_dna = True
        self.nsa_input2_is_dna = True
        self.nsa_input1_directionality = FIVEEND_TO_THIRDEND
        self.nsa_input2_directionality = FIVEEND_TO_THIRDEND
        self.icon_pixmap = QPixmap()

        self.setui()

        self.snct_list_of_first_base_label = [self.first_base_u_label, self.first_base_c_label, self.first_base_a_label, self.first_base_g_label]
        self.snct_list_of_second_base_label = [self.second_base_u_label, self.second_base_c_label, self.second_base_a_label, self.second_base_g_label]
        self.snct_list_of_third_base_label = [self.third_base_u_u_label, self.third_base_u_c_label, self.third_base_u_a_label, self.third_base_u_g_label,
                                              self.third_base_c_u_label, self.third_base_c_c_label, self.third_base_c_a_label, self.third_base_c_g_label,
                                              self.third_base_a_u_label, self.third_base_a_c_label, self.third_base_a_a_label, self.third_base_a_g_label,
                                              self.third_base_g_u_label, self.third_base_g_c_label, self.third_base_g_a_label, self.third_base_g_g_label]

        self.snct_dict_of_codon_label = {
            "UUU": self.uuu_label, "UCU": self.ucu_label, "UAU": self.uau_label, "UGU": self.ugu_label,
            "UUC": self.uuc_label, "UCC": self.ucc_label, "UAC": self.uac_label, "UGC": self.ugc_label,
            "UUA": self.uua_label, "UCA": self.uca_label, "UAA": self.uaa_label, "UGA": self.uga_label,
            "UUG": self.uug_label, "UCG": self.ucg_label, "UAG": self.uag_label, "UGG": self.ugg_label,
            "CUU": self.cuu_label, "CCU": self.ccu_label, "CAU": self.cau_label, "CGU": self.cgu_label,
            "CUC": self.cuc_label, "CCC": self.ccc_label, "CAC": self.cac_label, "CGC": self.cgc_label,
            "CUA": self.cua_label, "CCA": self.cca_label, "CAA": self.caa_label, "CGA": self.cga_label,
            "CUG": self.cug_label, "CCG": self.ccg_label, "CAG": self.cag_label, "CGG": self.cgg_label,
            "AUU": self.auu_label, "ACU": self.acu_label, "AAU": self.aau_label, "AGU": self.agu_label,
            "AUC": self.auc_label, "ACC": self.acc_label, "AAC": self.aac_label, "AGC": self.agc_label,
            "AUA": self.aua_label, "ACA": self.aca_label, "AAA": self.aaa_label, "AGA": self.aga_label,
            "AUG": self.aug_label, "ACG": self.acg_label, "AAG": self.aag_label, "AGG": self.agg_label,
            "GUU": self.guu_label, "GCU": self.gcu_label, "GAU": self.gau_label, "GGU": self.ggu_label,
            "GUC": self.guc_label, "GCC": self.gcc_label, "GAC": self.gac_label, "GGC": self.ggc_label,
            "GUA": self.gua_label, "GCA": self.gca_label, "GAA": self.gaa_label, "GGA": self.gga_label,
            "GUG": self.gug_label, "GCG": self.gcg_label, "GAG": self.gag_label, "GGG": self.ggg_label
        }

        self.snct_list_of_amino_acid = ["F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "O", "B", "H", "Q", "N", "K", "D", "E", "C", "J", "W", "R", "G"]

        self.snct_dict_of_amino_acid_label = {
            "F": [self.phenylalanine_label],
            "L": [self.leucine_label],
            "I": [self.isoleucine_label],
            "M": [self.methionine_label],
            "V": [self.valine_label],
            "S": [self.serine1_label, self.serine2_label],
            "P": [self.proline_label],
            "T": [self.threonine_label],
            "A": [self.alanine_label],
            "Y": [self.tyrosine_label],
            "O": [self.stop_ochre_label],
            "B": [self.stop_amber_label],
            "H": [self.histidine_label],
            "Q": [self.glutamine_label],
            "N": [self.asparagine_label],
            "K": [self.lysine_label],
            "D": [self.aspartic_acid_label],
            "E": [self.glutamic_acid_label],
            "C": [self.cysteine_label],
            "J": [self.stop_opal_label],
            "W": [self.tryptophan_label],
            "R": [self.arginine1_label, self.arginine2_label],
            "G": [self.glycine_label]
        }

        self.rnst_input_directionality = FIVEEND_TO_THIRDEND

        self.connect()

    def setui(self):
        self.setupUi(self)

        WindowTitle = f"Gene Expression Education Program v1.0.0"
        self.setWindowTitle(WindowTitle)
        window_icon = QIcon(ICON_DIR)
        self.setWindowIcon(window_icon)

        ICON_SIZE = 75
        self.icon_pixmap.load(ICON_DIR)
        self.icon_pixmap = self.icon_pixmap.scaled(QSize(ICON_SIZE, ICON_SIZE))
        self.nsa_icon_label.setPixmap(self.icon_pixmap)
        self.snct_icon_label.setPixmap(self.icon_pixmap)
        self.rnst_icon_label.setPixmap(self.icon_pixmap)

    def connect(self):
        """
        A method that connects pyqt5 signal
        :return:
        """
        # nsa
        self.nsa_input_lineEdit1.returnPressed.connect(lambda: self.nsa(LINE_EDIT_1))
        self.nsa_input_lineEdit2.returnPressed.connect(lambda: self.nsa(LINE_EDIT_2))
        self.nsa_input1DnaRadioButton.clicked.connect(lambda: self.nsa_is_dna_radio(LINE_EDIT_1, IS_DNA))
        self.nsa_input1RnaRadioButton.clicked.connect(lambda: self.nsa_is_dna_radio(LINE_EDIT_1, not IS_DNA))
        self.nsa_input2DnaRadioButton.clicked.connect(lambda: self.nsa_is_dna_radio(LINE_EDIT_2, IS_DNA))
        self.nsa_input2DnaRadioButton.clicked.connect(lambda: self.nsa_is_dna_radio(LINE_EDIT_1, not IS_DNA))
        self.nsa_input1_5_3RadioButton.clicked.connect(lambda: self.nsa_directionality_radio(LINE_EDIT_1, FIVEEND_TO_THIRDEND))
        self.nsa_input1_3_5RadioButton.clicked.connect(lambda: self.nsa_directionality_radio(LINE_EDIT_1, THIRDEND_TO_FIVEEND))
        self.nsa_input2_5_3RadioButton.clicked.connect(lambda: self.nsa_directionality_radio(LINE_EDIT_2, FIVEEND_TO_THIRDEND))
        self.nsa_input2_3_5RadioButton.clicked.connect(lambda: self.nsa_directionality_radio(LINE_EDIT_2, THIRDEND_TO_FIVEEND))

        # scnt
        self.u_checkBox_1.clicked.connect(self.snct)
        self.c_checkBox_1.clicked.connect(self.snct)
        self.a_checkBox_1.clicked.connect(self.snct)
        self.g_checkBox_1.clicked.connect(self.snct)
        self.u_checkBox_2.clicked.connect(self.snct)
        self.c_checkBox_2.clicked.connect(self.snct)
        self.a_checkBox_2.clicked.connect(self.snct)
        self.g_checkBox_2.clicked.connect(self.snct)
        self.u_checkBox_3.clicked.connect(self.snct)
        self.c_checkBox_3.clicked.connect(self.snct)
        self.a_checkBox_3.clicked.connect(self.snct)
        self.g_checkBox_3.clicked.connect(self.snct)

        # rnst
        self.rnst_input_lineEdit.returnPressed.connect(self.rnst)
        self.rnst_input_5_3RadioButton.clicked.connect(lambda: self.rnst_directionality_radio(FIVEEND_TO_THIRDEND))
        self.rnst_input_3_5RadioButton.clicked.connect(lambda: self.rnst_directionality_radio(THIRDEND_TO_FIVEEND))

    @pyqtSlot(int)
    def nsa(self, line_edit_num: int) -> None:
        """
        Nucleic Sequence Analyze method
        :param line_edit_num:
        :return:
        """
        if line_edit_num == LINE_EDIT_1:
            nucleic_sequence = self.nsa_input_lineEdit1.text()  # read lineEdit1
            dict_of_analyzed_data = analyze_nucleic_sequence(nucleic_sequence, self.nsa_input1_is_dna, self.nsa_input1_directionality)
            output = f"DNA1: {dict_of_analyzed_data['dna_strand1']}\n" \
                     f"RNA1: {dict_of_analyzed_data['rna_strand1']}\n" \
                     f"RNA2: {dict_of_analyzed_data['rna_strand2']}\n" \
                     f"DNA2: {dict_of_analyzed_data['dna_strand2']}\n" \
                     f"\n" \
                     f"RNA1: {dict_of_analyzed_data['rna_strand1']}\n" \
                     f"TNS1: {dict_of_analyzed_data['indented_target_area1']}\n" \
                     f"AAC1: {dict_of_analyzed_data['indented_translated_string1']}\n" \
                     f"\n" \
                     f"RNA2: {dict_of_analyzed_data['rna_strand2']}\n" \
                     f"TNS2: {dict_of_analyzed_data['indented_target_area2']}\n" \
                     f"AAC2: {dict_of_analyzed_data['indented_translated_string2']}"

            self.nsa_output1TextBrowser.setPlainText(output)  # print output

        elif line_edit_num == LINE_EDIT_2:
            nucleic_sequence = self.nsa_input_lineEdit2.text()  # read lineEdit2
            dict_of_analyzed_data = analyze_nucleic_sequence(nucleic_sequence, self.nsa_input2_is_dna, self.nsa_input2_directionality)
            output = f"DNA1: {dict_of_analyzed_data['dna_strand1']}\n" \
                     f"RNA1: {dict_of_analyzed_data['rna_strand1']}\n" \
                     f"RNA2: {dict_of_analyzed_data['rna_strand2']}\n" \
                     f"DNA2: {dict_of_analyzed_data['dna_strand2']}\n" \
                     f"\n" \
                     f"RNA1: {dict_of_analyzed_data['rna_strand1']}\n" \
                     f"TNS1: {dict_of_analyzed_data['indented_target_area1']}\n" \
                     f"AAC1: {dict_of_analyzed_data['indented_translated_string1']}\n" \
                     f"\n" \
                     f"RNA2: {dict_of_analyzed_data['rna_strand2']}\n" \
                     f"TNS2: {dict_of_analyzed_data['indented_target_area2']}\n" \
                     f"AAC2: {dict_of_analyzed_data['indented_translated_string2']}"

            self.nsa_output2TextBrowser.setPlainText(output)  # print output

    @pyqtSlot(int, bool)
    def nsa_is_dna_radio(self, line_edit_num: int, is_dna: bool) -> None:
        """
        A method that reads DNA/RNA radio button check
        :param line_edit_num:
        :param is_dna:
        :return:
        """
        if line_edit_num == LINE_EDIT_1:
            if is_dna:
                self.nsa_input1_is_dna = True
            else:
                self.nsa_input1_is_dna = False
        else:
            if is_dna:
                self.nsa_input2_is_dna = True
            else:
                self.nsa_input2_is_dna = False

    @pyqtSlot(int, str)
    def nsa_directionality_radio(self, line_edit_num: int, directionality: str) -> None:
        """
        A method that reads directionality radio button check
        :param line_edit_num:
        :param directionality:
        :return:
        """
        if directionality in [FIVEEND_TO_THIRDEND, THIRDEND_TO_FIVEEND]:
            if line_edit_num == LINE_EDIT_1:
                self.nsa_input1_directionality = directionality

            elif line_edit_num == LINE_EDIT_2:
                self.nsa_input2_directionality = directionality

    @pyqtSlot()
    def snct(self) -> None:
        """
        Specific Nucleo base Combine and Translation method
        :return:
        """

        LIST_OF_BASE = ['U', 'C', 'A', 'G']
        LIST_OF_CODON = ["UUU", "UCU", "UAU", "UGU",
                         "UUC", "UCC", "UAC", "UGC",
                         "UUA", "UCA", "UAA", "UGA",
                         "UUG", "UCG", "UAG", "UGG",
                         "CUU", "CCU", "CAU", "CGU",
                         "CUC", "CCC", "CAC", "CGC",
                         "CUA", "CCA", "CAA", "CGA",
                         "CUG", "CCG", "CAG", "CGG",
                         "AUU", "ACU", "AAU", "AGU",
                         "AUC", "ACC", "AAC", "AGC",
                         "AUA", "ACA", "AAA", "AGA",
                         "AUG", "ACG", "AAG", "AGG",
                         "GUU", "GCU", "GAU", "GGU",
                         "GUC", "GCC", "GAC", "GGC",
                         "GUA", "GCA", "GAA", "GGA",
                         "GUG", "GCG", "GAG", "GGG"]

        list_of_first_base = []
        list_of_second_base = []
        list_of_third_base = []
        list_of_first_base_boolean = [self.u_checkBox_1.isChecked(), self.c_checkBox_1.isChecked(), self.a_checkBox_1.isChecked(), self.g_checkBox_1.isChecked()]
        list_of_second_base_boolean = [self.u_checkBox_2.isChecked(), self.c_checkBox_2.isChecked(), self.a_checkBox_2.isChecked(), self.g_checkBox_2.isChecked()]
        list_of_third_base_boolean = [self.u_checkBox_3.isChecked(), self.c_checkBox_3.isChecked(), self.a_checkBox_3.isChecked(), self.g_checkBox_3.isChecked()]

        # check selected base
        for index, base in enumerate(LIST_OF_BASE, start=1):
            if list_of_first_base_boolean[index]:
                list_of_first_base.append(base)

            if list_of_second_base_boolean[index]:
                list_of_second_base.append(base)

            if list_of_third_base_boolean[index]:
                list_of_third_base.append(base)

        list_of_selected_codon, list_of_amino_acid = match_digit_base_to_codon_and_amino_acid(list_of_first_base_boolean, list_of_second_base_boolean, list_of_third_base_boolean)

        # reset base label color
        for label in self.snct_list_of_first_base_label:
            label.setStyleSheet("background-color: rgba(255, 0, 0, 0); border: 1px solid black;")
            # label.setStyleSheet("background:transparent")

        for label in self.snct_list_of_second_base_label:
            label.setStyleSheet("background-color: rgba(255, 0, 0, 0); border: 1px solid black;")
            # label.setStyleSheet("background:transparent")

        for label in self.snct_list_of_third_base_label:
            label.setStyleSheet("background-color: rgba(255, 0, 0, 0); border: 1px solid black;")

        # reset codon label color
        for codon in LIST_OF_CODON:
            label = self.snct_dict_of_codon_label[codon]
            label.setStyleSheet("background-color: rgba(255, 0, 0, 0); border: 1px solid black;")

        # reset amino acid label color
        for amino_acid in self.snct_list_of_amino_acid:
            for label in self.snct_dict_of_amino_acid_label[amino_acid]:
                label.setStyleSheet("background-color: rgba(255, 0, 0, 0); border: 1px solid black;")

        # highlight bases label
        for first_base in list_of_first_base:
            self.snct_list_of_first_base_label[LIST_OF_BASE.index(first_base)].setStyleSheet("background-color: lightgreen; border: 1px solid black;")

        for second_base in list_of_second_base:
            self.snct_list_of_second_base_label[LIST_OF_BASE.index(second_base)].setStyleSheet("background-color: lightgreen; border: 1px solid black;")

        for third_base in list_of_third_base:
            self.snct_list_of_third_base_label[LIST_OF_BASE.index(third_base) + 4 * 0].setStyleSheet("background-color: lightgreen; border: 1px solid black;")
            self.snct_list_of_third_base_label[LIST_OF_BASE.index(third_base) + 4 * 1].setStyleSheet("background-color: lightgreen; border: 1px solid black;")
            self.snct_list_of_third_base_label[LIST_OF_BASE.index(third_base) + 4 * 2].setStyleSheet("background-color: lightgreen; border: 1px solid black;")
            self.snct_list_of_third_base_label[LIST_OF_BASE.index(third_base) + 4 * 3].setStyleSheet("background-color: lightgreen; border: 1px solid black;")

        # highlight codon label
        for codon_ in list_of_selected_codon:
            self.snct_dict_of_codon_label[codon_].setStyleSheet("background-color: lightgreen; border: 1px solid black;")

        # highlight amino acid label
        for amino_acid in list_of_amino_acid:
            list_of_amino_acid_label = self.snct_dict_of_amino_acid_label[amino_acid]
            for label in list_of_amino_acid_label:
                label.setStyleSheet("background-color: yellow; border: 1px solid black;")

    @pyqtSlot(str)
    def rnst_directionality_radio(self, directionality: str) -> None:
        """
        A method that reads directionality radio button check
        :param directionality:
        :return:
        """
        if directionality in [FIVEEND_TO_THIRDEND, THIRDEND_TO_FIVEEND]:
            self.rnst_input_directionality = directionality

    @pyqtSlot()
    def rnst(self) -> None:
        """
        Repeated Nucleic Sequence Translation method
        :return:
        """
        nucleic_sequence = self.rnst_input_lineEdit.text()
        list_of_analysis_pair = combine_all_amino_acids_from_repeated_nucleic_sequence(nucleic_sequence, self.rnst_input_directionality)
        output = ""
        column = 0
        for analysis_pair in list_of_analysis_pair:
            codon, amino_acid = analysis_pair
            output = f"{output}{codon}: {amino_acid}  "
            column += 1

            if column == 8:
                output = f"{output}\n"
                column %= 8

        self.rnst_outputTextBrowser.setPlainText(output)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    myApp = MainWindow()
    myApp.show()
    sys.exit(app.exec_())  # run GUI

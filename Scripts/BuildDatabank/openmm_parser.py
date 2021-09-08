from xml.etree import ElementTree as ET

def openmm_input(input, type):

    if type == 'inp':

        return read_inp_input(input)

    elif type == 'xml':

        return read_xml_input(input)


def read_inp_input(input):

    try:

        filename = open(input, 'r')
        print(filename.name)

    except FileNotFoundError:

        print('File does not exist')

    with open(filename.name, 'r') as f:

        tmp_type = []
        tmp_value = []


        for line in f:

            split = line.split()

            if len(split) > 2:

                tmp_type.append(split[0])
                tmp_value.append(split[2])

    return [tmp_type, tmp_value]

def read_xml_input(input):

    try:

        filename = open(input, 'r')

    except FileNotFoundError:

        print('File does not exist')

    tree = ET.parse(filename.name)
    temp = [value for key, value in tree.getroot()[1].attrib.items() if 'temperature' in key.lower()]

    return(['temp'],[temp])






class openmmParser:

    def __init__(self, filename, type):
        [properties, values] = openmm_input(filename, type)
        print(properties)
        self.temperature = float(values[properties.index('temp')])
 #       self.pressure = values[properties.index('p_ref')]
 #       self.ptype = values[properties.index('p_type')]
 #       self.dt = values[properties.index('dt')]

    def setTempurature(self, temp):

        self.temperature = temp

    def setPressure(self, pressure):

        self.pressure = pressure

    def setPType(self, ptype):

        self.ptype = ptype

    def setDt(self, dt):

        self.dt = dt








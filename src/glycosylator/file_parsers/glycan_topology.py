import copy
import sqlite3

# TODO: WIP


def read_glycan_topology(file_path):
    """read glycan.top file"""
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines()]
        lines = [line.split("!")[0] for line in lines]
        lines = [line for line in lines if line != ""]
    residue = {}
    connect_topology = {}
    nbr_units = 0
    for line in lines:
        field, *value = line.split()
        if field == "RESI":
            if residue:
                residue["#UNIT"] = nbr_units
                connect_topology[resname] = copy.copy(residue)
            residue["UNIT"] = []
            resname = line[1]
            nbr_units = 0
        elif field == "UNIT":
            res_name, *path = value

            read_unit(line.split(), residue)
            nbr_units += 1

            # if len(unit) > 2:

    #     residue["UNIT"].append([unit[1], unit[2], unit[3:]])
    # else:
    #     residue["UNIT"].append([unit[1], " ", []])

    residue["#UNIT"] = nbr_units
    connect_topology[resname] = copy.copy(residue)


def import_connectivity_topology(filename):
    """Import connectivity topology from sql database
    This function will initialize connect_topology
    Parameters:
        filename: path to database
    """
    try:
        conn = sqlite3.connect(filename)
    except:
        print(f"Error while connecting to the database {filename}")
        return -1
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM glycans")
    glycans = cursor.fetchall()

    connect_topology = {}
    for glycan in glycans:
        name, tree = glycan
        residue = {}
        residue["UNIT"] = []
        nbr_unit = 0
        for unit in tree.split("|"):
            unit = unit.split(" ")
            nbr_unit += 1
            if len(unit) > 2:
                residue["UNIT"].append([unit[0], unit[1], unit[2:]])
            else:
                residue["UNIT"].append([unit[0], " ", []])

        residue["#UNIT"] = nbr_unit
        connect_topology[name] = residue
    return connect_topology


def export_connectivity_topology(filename, connect_topology):
    """Export connectivity topology to sql database
    This function will
    """
    try:
        conn = sqlite3.connect(filename)
    except:
        print(f"Error while connecting to the database {filename}")
        return -1
    cursor = conn.cursor()
    tn = "glycans"
    gn = "glycan_name"
    gt = "glycan_tree"
    cursor.execute(f"DROP TABLE IF EXISTS {tn}")
    cursor.execute(f"CREATE TABLE {tn} ({gn} text, {gt} text)")

    for key in connect_topology.keys():
        units = connect_topology[key]["UNIT"]
        glycan = []
        for unit in units:
            v = []
            v.extend(unit[0:2])
            v.extend(unit[2])
            glycan.append(" ".join(v))
        glycan = "|".join(glycan)

        cursor.execute(f"INSERT INTO {tn} VALUES ('{key}', '{glycan}')")

    conn.commit()
    conn.close()


def add_glycan_to_connectivity_topology(self, name, linkage_paths, overwrite=True):
    """Add new glycan to connect_topology dictionary
    Parameters:
        name: name of new glycan
        connect_tree: dictionary with connectivity tree
        overwrite: should an existing glycan be overwritten
    """
    if name in self.connect_topology and not overwrite:
        print(
            "Glycan with same name "
            + name
            + "already exists. Please change name or allow overwritting"
        )
        return -1
    self.connect_topology[name] = linkage_paths


def read_unit(unit, residue):
    if len(unit) > 2:
        residue["UNIT"].append([unit[1], unit[2], unit[3:]])
    else:
        residue["UNIT"].append([unit[1], " ", []])

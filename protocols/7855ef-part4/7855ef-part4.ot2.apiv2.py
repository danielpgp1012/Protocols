import math
from opentrons import protocol_api
import os
import csv

metadata = {
    'protocolName': 'Agriseq Library Prep Part 4 - Pooling',
    'author': 'Rami Farawi <rami.farawi@opentrons.com>',
    'source': 'Custom Protocol Request',
    'apiLevel': '2.7'
}


def run(protocol):

    [num_samp, p20_mount, p300_mount,
        reset_tipracks] = get_values(  # noqa: F821
        "num_samp", "p20_mount", "p300_mount", "reset_tipracks")

    if not 1 <= num_samp <= 288:
        raise Exception("Enter a sample number between 1-288")

    # Tip tracking between runs
    if not protocol.is_simulating():
        file_path = '/data/csv/tiptracking.csv'
        file_dir = os.path.dirname(file_path)
        # check for file directory
        if not os.path.exists(file_dir):
            os.makedirs(file_dir)
        # check for file; if not there, create initial tip count tracking
        if not os.path.isfile(file_path):
            with open(file_path, 'w') as outfile:
                outfile.write("0, 0\n")

    tip_count_list = []
    if protocol.is_simulating():
        tip_count_list = [0, 0]
    elif reset_tipracks:
        tip_count_list = [0, 0]
    else:
        with open(file_path) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            tip_count_list = next(csv_reader)

    num_one = int(tip_count_list[0])
    num_two = int(tip_count_list[1])

    # load labware
    pool_plate1 = protocol.load_labware('customendura_96_wellplate_200ul', '1',
                                        label='Pool Plate 1')
    pool_plate2 = protocol.load_labware('customendura_96_wellplate_200ul', '2',
                                        label='Pool Plate 2')
    reaction_plates = [protocol.load_labware('customendura_96_wellplate_200ul',
                       str(slot), label='Reaction Plate')
                       for slot in [4, 5, 6]]
    tiprack20 = [protocol.load_labware('opentrons_96_filtertiprack_20ul',
                                       str(slot))
                 for slot in [9, 10, 11]]

    tiprack200 = [protocol.load_labware(
                    'opentrons_96_filtertiprack_200ul', '8')]

    # load instruments
    p20 = protocol.load_instrument('p20_single_gen2', p20_mount,
                                   tip_racks=tiprack20)
    p300 = protocol.load_instrument('p300_single_gen2', p300_mount,
                                    tip_racks=tiprack200)

    def pickup(pip):
        try:
            pip.pick_up_tip()
        except protocol_api.labware.OutOfTipsError:
            protocol.pause("Replace all 20 ul tip racks on Slots 9, 10 and 11")
            pip.reset_tipracks()
            pip.pick_up_tip()

    # pool each row of plates
    num_col = math.ceil(num_samp/8)
    pool_counter = 0
    vol_counter = 0
    num_plates = math.ceil(num_samp/96)
    for i, plate in enumerate(reaction_plates[:num_plates]):
        if i == 0:
            length_row = num_col
        elif i == 1:
            length_row = num_col - 12
        else:
            length_row = num_col - 24
        wells_by_row = [well for row in reaction_plates[i].rows()
                        for well in row[:length_row]]
        for well in wells_by_row:
            pickup(p20)
            p20.aspirate(5, well)
            p20.dispense(5, pool_plate1.wells()[pool_counter])
            p20.blow_out()
            vol_counter += 5
            if vol_counter % 60 == 0:
                pool_counter += 1
                protocol.comment('\nNEXT POOL WELL\n')
            p20.return_tip()
            protocol.comment('\n')

    for s, d in zip(pool_plate1.wells()[:pool_counter], pool_plate2.wells()):
        pickup(p300)
        p300.mix(2, 60)
        p300.transfer(45, s, d, new_tip='never')
        p300.return_tip()

    # write updated tipcount to CSV
    new_tip_count = str(num_one)+", "+str(num_two)+"\n"
    if not protocol.is_simulating():
        with open(file_path, 'w') as outfile:
            outfile.write(new_tip_count)

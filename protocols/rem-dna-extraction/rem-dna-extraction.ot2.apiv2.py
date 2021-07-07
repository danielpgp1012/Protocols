# initial disposition
# slot1:magnetic module (don't forget to turn it on from its back)
# slot2:reservoir A2-3 ethanol(full for full plate) A5 EB(2/3rd is enough) -->this protocol adjusts the height of pipetting according to the numb of columns.
# slot3:tiprack1  WILL BE REUSED. PUT ONLY THE NUMBER OF COLUMNS THAT YOU NEED, STARTING FROM LEFT
# slot4:tiprack2  WILL BE REUSED. PUT ONLY THE NUMBER OF COLUMNS THAT YOU NEED, STARTING FROM LEFT
# slot5:CHANGES: 1st: plate with samples, so that magnet doesn't affect binding to beads;
# 2nd: plate to trash samples;
# AND 3rd plate for DNA samples
# slot6:tiprack3  *WILL BE REUSED* if less than 12 columns, the rest of the tiprack should be filled because it will be used for next steps
# slot7:tiprack4 #full tiprack
# slot8:tiprack5 #full tiprack

# slot9:tiprack6   should not be required (empty)
# slot10:tiprack7  should not be required (empty)
# slot11:tiprack8  should not be required (empty)

# to be adjusted BEFORE LOADING TO OPENTRONS:
# line 57: number of columns -->VERY important
# line 60: z_offset depending on sample viscosity -->mostly is important for removing SN with lysis buffer and avoid taking up the beads
# line 64: volume of starting sample
# line 65: volume of beads (for total volume calculation)

# Number of Pauses that have to be manually resumed: 4
# At start
# When putting plate on magnet module
# When drying plate before adding EB.
# When putting plate at 60°C for DNA elution from beads

# "protocol.comment(print(used_tipracks, "tipracks have been used.")" comments are for the debugging process
from opentrons import protocol_api
from opentrons.types import Point
import random

metadata = {
    "apiLevel": "2.8",
    "protocolName": "DNA_extraction_beads",
    "author": "REM",
    "description": "DNA extraction 96w plate",
}


def run(protocol: protocol_api.ProtocolContext):
    # columns
    column_no = 4
    # Number of mixing times for lysate + beads
    lysate_mix_reps = 25
    # lysate and bead volumes
    [v_sample, v_beads] = [50, 20]
    # ethanol volume in uL
    v_etoh = 120
    # final elution buffer
    v_EB = 50
    # binding time in minutes
    bead_binding_time = 10
    mag_binding_time = 5
    bead_detaching_time_in_s = 5
    # beads
    bead_tube = protocol.load_labware(
        "opentrons_24_tuberack_eppendorf\
        _1.5ml_safelock_snapcap",
        8,
    )["A1"]
    # ethanol
    reagents = protocol.load_labware("nest_12_reservoir_15ml", 2)
    # tip racks
    n_tiprack = 3

    # multi-channel tiprack
    slots = [3, 4, 6, 7]
    assert n_tiprack <= 4
    m_tipacks = [
        protocol.load_labware("opentrons_96_tiprack_300ul", i)
        for i in [3, 4, 6, 7, 8][:n_tiprack]
    ]

    # single-channel tiprack
    s_tiprack = protocol.load_labware("opentrons_96_tiprack_300ul", 9)

    # reaction plate (single one that will be moved by user)
    pre_l_plate = protocol.load_labware("96w_pcr_plate2", 5)
    mag_mod: protocol_api.MagneticModuleContext = protocol.load_module(
        "magnetic module", 1
    )
    r_plate = mag_mod.load_labware("96w_pcr_plate2")

    # pipettes
    # multi-channel
    multi_300: protocol_api.InstrumentContext = protocol.load_instrument(
        "p300_multi_gen2", mount="right", tip_racks=m_tipacks
    )
    tip300_count = 0
    tip300_max = 12 * n_tiprack

    def multi_tip_pick_up():
        """Refillable tip pickup.
        If no tips left, protocol will
        pause and ask for refill."""
        nonlocal tip300_count
        if tip300_count == tip300_max:
            protocol.pause(
                f"Refill 300ul tipracks in slots {slots} before resuming."
            )
            multi_300.reset_tipracks()
            tip300_count = 0
        tip300_count += 1
        multi_300.pick_up_tip()

    # single channel
    single_300: protocol_api.InstrumentContext = protocol.load_labware(
        "p300_single_gen2", mount="left", tip_racks=[s_tiprack]
    )

    # lights on
    protocol.set_rail_lights(False)
    protocol.set_rail_lights(True)

    # step 1 magnet
    mag_mod.disengage()
    protocol.pause(
        "Put plate on slot 5 (not magnet) and start. It should \
            contain up to 100uL of sample and 20uL of magnetic beads."
    )
    # step 2 add beads
    single_300.flow_rate.aspirate = 67
    single_300.flow_rate.dispense = 120

    for i in range(column_no):
        single_300.distribute(
            v_beads,
            bead_tube,
            pre_l_plate.columns()[i],
            air_gap=20,
            new_tip="always",
            mix_before=(4, 300),
        )

    # step 3 mix with multi-channel
    multi_300.flow_rate.aspirate = 67
    multi_300.flow_rate.dispense = 120

    # 2 ok for 50uL+20uL, 4 for 120uL total
    multi_300.well_bottom_clearance.aspirate = 3
    multi_300.well_bottom_clearance.dispense = 2

    # mix and return tips
    for i in range(column_no):
        multi_300.pick_up_tip()
        multi_300.mix(
            lysate_mix_reps,
            min((v_beads + v_sample) * 0.7, multi_300.max_volume),
            pre_l_plate.columns()[i][0],
        )
        multi_300.blow_out()
        multi_300.return_tip()
    multi_300.reset_tipracks()

    # step 4 bead binding
    protocol.delay(
        minutes=bead_binding_time,
        msg=f"RT: {str(bead_binding_time)} minutes for\
             DNA binding to beads.",
    )
    protocol.pause(
        "Place sample plate on magnetic module.\
        Place a trash plate on slot 5"
    )

    # step 5 engage magnet
    mag_mod.engage(height=15)
    protocol_api.delay(
        minutes=mag_binding_time,
        msg=f"Wait {str(mag_binding_time)} minutes for\
             beads to bind to magnet.Make sure beads\
                  are visibly sticking to the wall!",
    )

    # step 6 transfer EDM supernatant to waste
    v_asp_sn = (v_sample + v_beads) * 1.4
    multi_300.flow_rate.aspirate = 30
    multi_300.flow_rate.dispense = 300
    multi_300.well_bottom_clearance.dispense = 5
    multi_300.well_bottom_clearance.aspirate = 1.7
    x_left = -1.1
    x_right = 1.1
    for i in range(column_no):
        multi_tip_pick_up()
        # offset x location dependening on magnet location
        not_bead_loc = x_right if (i + 1) % 2 == 0 else x_left
        multi_300.transfer(
            v_asp_sn,
            r_plate.columns()[i].move(Point(x=not_bead_loc)),
            pre_l_plate.columns()[i],
            air_gap=20,
        )
        multi_300.drop_tip()
    # step 7 disengage magnet
    mag_mod.disengage()
    multi_300.flow_rate.aspirate = 200
    multi_300.flow_rate.dispense = 300
    multi_300.well_bottom_clearance.dispense = 6

    # step 8: ethanol wash. step 9: mix
    def mix_at_custom_speed_at_random_offsets(
        reps,  # repetitions
        vol,  # volumes
        dest,  # destination point
        asp_r,  # aspirationrate
        disp_r=None,  # dispense rate defaults to equal asp_r
        x_limits=(-1.1, 1.1),  # min-max pipette locations
        y_limits=(-0.8, 0.8),
        z_limits=(-1, 3),
    ):
        """Wrapper of mixing functionality.
        We customize aspiration and dispense speeds.
        reps: repetitions
        vol: volume of mixing
        dest: destination of mixing
        asp_r,disp_r: multiplication factor of base flow rate. If 4, then flow rate will be 4 times base flow rate
        x_limits, y_limits, z_limits: range of locations for pipette within well"""
        if not disp_r:
            disp_r = asp_r
        for i in range(reps):
            x = round(random.uniform(x_limits[0], x_limits[1]), 1)
            y = round(random.uniform(y_limits[0], y_limits[1]), 1)
            z = round(random.uniform(z_limits[0], z_limits[1]), 1)
            offset = Point(x, y, z)
            multi_300.aspirate(vol, dest.move(offset), asp_r, disp_r)
            multi_300.dispense(vol, dest.move(offset), asp_r, disp_r)

    for i in range(column_no):
        multi_tip_pick_up()
        # transfer ethanol to reaction plate
        multi_300.transfer(
            v_etoh, reagents["A2"], r_plate.columns()[i], rate=0.5, air_gap=30
        )
        bead_location = x_right if (i + 1) % 2 == 0 else x_right
        dest = r_plate.columns()[i]
        # mix at bead pellet
        mix_at_custom_speed_at_random_offsets(
            5,  # repetitions
            v_etoh,  # volume
            dest,  # destination
            5,  # times base flow rate
            x_limits=(bead_location, bead_location),  # focus on bead pellet
        )
        # generalized mix
        mix_at_custom_speed_at_random_offsets(10, v_etoh, dest, 4)
        protocol.delay(seconds=2)
        multi_300.blow_out(dest.top(z=-1))
        multi_300.return_tip()
    # step 10 engage magnet
    mag_mod.engage(height=12)

    # step 11
    protocol.delay(
        seconds=bead_detaching_time_in_s,
        msg="Turn on heatblock at 60°C.\
             Wait until magnet acts properly.\
             make sure it is sticking to wall, or pause manually. If needed remove tips from trash.",
    )

import math
import matplotlib.pyplot as plt

start_time = 0
end_time_hr = 40
calculation_time_step_hr = 1/3600

#time_steps = list(range(start_time, end_time_hr, calculation_time_step_hr))
# take the time steps list for the calculation
sequence = list(range(int(start_time / calculation_time_step_hr), int(end_time_hr / calculation_time_step_hr) + 1))
time_steps = [x * calculation_time_step_hr for x in sequence]

# Question input parameters
Qb = 2
Qs = 45
Ta = 9
nf = 12
Qo_max = 12.984953539
loss_coefficient = 0.55
downstream_wl = 324

# Quadratic equation coefficients as given in the question
a = 1450
b = 5645

# create a list for Qin and Qout result store
Qin_results = []
Qout_results = []
reservoir_wl_results = []
reservoir_wl = 0
effective_WH = reservoir_wl - downstream_wl
cumulative_storage = 0
cumulative_storage_results = []


# Range of diameters to check
min_diameter = 1.259
max_diameter = 2.0
step_diameter = 0.001

best_diameter = None
best_max_discharge = float('-inf')


for assumed_dia in range(int(min_diameter / step_diameter), int(max_diameter / step_diameter) + 1):
    assumed_dia = assumed_dia * step_diameter
    area_calc = math.pi * (assumed_dia ** 2) / 4

    # create a list for Qin and Qout result store
    Qin_results = []
    Qout_results = []
    reservoir_wl_results = []
    reservoir_wl = 0
    effective_WH = reservoir_wl - downstream_wl
    cumulative_storage = 0
    cumulative_storage_results = []

    for t in time_steps:

        # inflow calculation
        Qin = Qb + (Qs - Qb) * ((t / Ta) * math.exp(1 - (t / Ta))) ** nf
        Qin_results.append(Qin)

        # Outflow calculation
        if effective_WH <= 0:
            Qout = 0

        else:
            Qout = loss_coefficient * area_calc * math.sqrt(2 * 9.81 * effective_WH)

        Qout_results.append(Qout)

        # net inflow calculation
        net_inflow = Qin - Qout

        # store volume calculation
        stored_volume = net_inflow * 3600 * calculation_time_step_hr

        # update of cumulative storage
        cumulative_storage += stored_volume
        cumulative_storage_results.append(cumulative_storage)

        # quadratic equation last coefficient (change with time)
        c = -cumulative_storage

        # calculation of the reservoir level
        h1 = ((-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)) + 322
        h2 = ((-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)) + 322

        reservoir_wl = max(h1, h2)
        water_height = reservoir_wl - 322
        effective_WH = reservoir_wl - downstream_wl

        reservoir_wl_results.append(reservoir_wl)

    if max(Qout_results) >= Qo_max:
        best_max_discharge = max(Qout_results)
        best_diameter = assumed_dia - step_diameter
        break


# run the code for selected diameter.

area_calc_final = math.pi * (best_diameter ** 2) / 4

# create a list for Qin and Qout result store
Qin_results = []
Qout_results = []
reservoir_wl_results = []
reservoir_wl = 0
effective_WH = reservoir_wl - downstream_wl
water_height_results = []
cumulative_storage = 0
cumulative_storage_results = []

for t in time_steps:

    # inflow calculation
    Qin = Qb + (Qs - Qb) * ((t / Ta) * math.exp(1 - (t / Ta))) ** nf
    Qin_results.append(Qin)

    # Outflow calculation
    if effective_WH <= 0:
        Qout = 0

    else:
        Qout = loss_coefficient * area_calc_final * math.sqrt(2 * 9.81 * effective_WH)

    Qout_results.append(Qout)

    # net inflow calculation
    net_inflow = Qin - Qout

    # store volume calculation
    stored_volume = net_inflow * 3600 * calculation_time_step_hr

    # update of cumulative storage
    cumulative_storage += stored_volume
    cumulative_storage_results.append(cumulative_storage)

    # quadratic equation last coefficient (change with time)
    c = -cumulative_storage

    # calculation of the reservoir level
    h1 = ((-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)) + 322
    h2 = ((-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)) + 322

    reservoir_wl = max(h1, h2)
    water_height = reservoir_wl - 322
    effective_WH = reservoir_wl - downstream_wl

    water_height_results.append(water_height)

    reservoir_wl_results.append(reservoir_wl)

# Output results
print("Best Diameter:", best_diameter)
print("Best Maximum Discharge:", best_max_discharge)
print("Maximum water level: ",max(reservoir_wl_results), "m asl")
print("Maximum water height / stage : ",max(water_height_results), "m asl")
print("Maximum water level reach at:",(reservoir_wl_results.index(max(reservoir_wl_results)))*calculation_time_step_hr, " hr")
print("Storage Capacity: ", max(cumulative_storage_results))
print("Maximum dischrge obtained: ", max(Qout_results))


# Create the first subplot for Qin and Qout
fig, ax1 = plt.subplots()

ax1.plot(time_steps, Qin_results, label='Inflow (m3/s) ', color='blue')
ax1.plot(time_steps, Qout_results, label='Discharge (m3/s)', color='green')
ax1.set_xlabel('Time(hr)')
ax1.set_ylabel('Inflow and Outflow (m3/s)', color='black')
ax1.tick_params('y', colors='black')
ax1.legend(loc='upper left')

# Create the second subplot for Reservoir WL with a different y-axis
ax2 = ax1.twinx()
ax2.plot(time_steps, water_height_results, label='Reservoir Water Height / Stage (m)', color='red')
ax2.set_ylabel('Reservoir Water Height / Stage (m)', color='black')
ax2.tick_params('y', colors='black')
ax2.legend(loc='upper right')

# Mark the maximum water level with a horizontal dotted line and label
max_water_level = max(water_height_results)
ax2.axhline(y=max_water_level, color='black', linestyle='--', label='Max Water Level')
ax2.text(end_time_hr / 2, max_water_level, f'Maximum Water Level / Stage: {max_water_level:.2f} m', color='black')

# Mark the time for the maximum water level with a vertical dotted line and label
time_max_water_level = water_height_results.index(max_water_level) * calculation_time_step_hr
ax2.axvline(x=time_max_water_level, color='orange', linestyle='--', label='Time for Max Water Level / Stage')
#ax2.text(time_max_water_level, max_water_level, f'Time for Max Water Level: {time_max_water_level:.2f} hr', color='black', rotation=90, va='bottom')

# Mark the maximum discharge with a horizontal dotted line and label
max_discharge = max(Qout_results)
ax1.axhline(y=max_discharge, color='purple', linestyle='--', label='Max Discharge')
ax1.text(end_time_hr / 2, max_discharge, f'Maximum Discharge: {max_discharge:.2f} m3/s', color='purple')

# Show the legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left')

plt.show()

# Create the first figure for Qin
plt.figure()
plt.plot(time_steps, Qin_results, label='Inflow (m3/s) ', color='blue')
plt.axhline(y=max(Qin_results), color='black', linestyle='--', label='Max Inflow')
plt.text(end_time_hr / 2, max(Qin_results), f'Max Inflow: {max(Qin_results):.2f} m3/s', color='black')
plt.axvline(x=time_steps[Qin_results.index(max(Qin_results))], color='orange', linestyle='--', label='Time for Max Inflow')
plt.text(time_steps[Qin_results.index(max(Qin_results))], 0, f'Time to Max Inflow: {time_steps[Qin_results.index(max(Qin_results))]:.2f} hr', color='black', rotation=90, va='bottom')
plt.xlabel('Time(hr)')
plt.ylabel('Inflow (m3/s)')
plt.title('Inflow (Qin)')
plt.legend()
plt.show()

# Create the second figure for Qout
plt.figure()
plt.plot(time_steps, Qout_results, label='Discharge (m3/s)', color='green')
plt.axhline(y=max(Qout_results), color='black', linestyle='--', label='Max Discharge')
plt.text(end_time_hr / 2, max(Qout_results), f'Max Discharge: {max(Qout_results):.2f} m3/s', color='black')
plt.axvline(x=time_steps[Qout_results.index(max(Qout_results))], color='orange', linestyle='--', label='Time for Max Discharge')
plt.text(time_steps[Qout_results.index(max(Qout_results))], 0, f'Time to Max Discharge: {time_steps[Qout_results.index(max(Qout_results))]:.2f} hr', color='black', rotation=90, va='bottom')
plt.xlabel('Time(hr)')
plt.ylabel('Discharge (m3/s)')
plt.title('Discharge (Qout)')
plt.legend()
plt.show()

# Create the third figure for Reservoir Water Height
plt.figure()
plt.plot(time_steps, water_height_results, label='Reservoir Water Level ASL (m)', color='red')
plt.axhline(y=max(water_height_results), color='black', linestyle='--', label='Max Water Level')
plt.text(end_time_hr / 2, max(water_height_results), f'Max Water Height / Stage: {max(water_height_results):.2f} m', color='black')
plt.axvline(x=time_steps[water_height_results.index(max(water_height_results))], color='orange', linestyle='--', label='Time for Max Water Level')
plt.text(time_steps[water_height_results.index(max(water_height_results))], 0, f'Time to Max Water Level: {time_steps[water_height_results.index(max(water_height_results))]:.2f} hr', color='black', rotation=90, va='bottom')
plt.xlabel('Time(hr)')
plt.ylabel('Reservoir Water Height / Stage (m)')
plt.title('Reservoir Water Height / Stage')
plt.legend()
plt.show()


# Create the third figure for Reservoir Water Level
plt.figure()
plt.plot(time_steps, reservoir_wl_results, label='Reservoir Water Level ASL (m)', color='brown')
plt.axhline(y=max(reservoir_wl_results), color='black', linestyle='--', label='Max Water Level')
plt.text(end_time_hr / 2, max(reservoir_wl_results), f'Max Water Level: {max(reservoir_wl_results):.2f} m', color='black')
plt.axvline(x=time_steps[reservoir_wl_results.index(max(reservoir_wl_results))], color='black', linestyle='--', label='Time for Max Water Level')
plt.text(time_steps[reservoir_wl_results.index(max(reservoir_wl_results))], 0, f'Time to Max Water Level: {time_steps[reservoir_wl_results.index(max(reservoir_wl_results))]:.2f} hr', color='black', rotation=90, va='bottom')
plt.xlabel('Time(hr)')
plt.ylabel('Reservoir Water Level (m)')
plt.title('Reservoir Water Level')
plt.legend()
plt.show()
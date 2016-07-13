import shelve
import matplotlib.pyplot as plt

__author__ = 'Kyle Vitautas Lopin'

result_shelf = shelve.open("Results1")

results_eig = result_shelf["results eig"]
results_svd = result_shelf["results svd"]
result_shelf.close()

states = []
labels = []
for i in range(9):
    states.append([])
    labels.append("line " + str(i))

print dir(results_eig[0])
print results_eig[0]
i2 = 0

for result in results_eig:
    # print result.voltage
    # print result.steady_state[5], result.steady_state[8]
    # print 'i2: ', i2, result.voltage
    i2 += 1
    # print result.steady_state
    for i, state in enumerate(result.steady_state):

        states[i].append(state)
# print states

# Na_transport_figure = plt.figure()
# Na_transport_axis = Na_transport_figure.add_subplot(111)
# Na_transport_axis.plot(voltages, Na_transport0, voltages, Na_transport1, voltages, Na_transport2)
voltages = range(-150, 100, 10)
state_fig = plt.figure()
state_axis = state_fig.add_subplot(111)
#state_axis.plot(voltages, states[0])
for state in states:
    state_axis.plot(voltages, state)
    print state


state_axis.legend(labels)
state_axis.set_yscale('log')
plt.show()
# plt.legend(handles=range(0, 9))
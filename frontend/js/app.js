// Constants
const NUM_REGISTERS = 8;
let registers = Array(NUM_REGISTERS).fill(0);
let memory = [
    0, 10,
    1, 20,
    2, 30
]; // example memory

// Update the tables
function updateTables() {
    const regTable = document.getElementById("registerTable");
    regTable.innerHTML = "<tr><th>Register</th><th>Value</th></tr>";
    for (let i = 0; i < NUM_REGISTERS; i++) {
        regTable.innerHTML += `<tr><td>R${i}</td><td>${registers[i]}</td></tr>`;
    }

    const memTable = document.getElementById("memoryTable");
    memTable.innerHTML = "<tr><th>Address</th><th>Value</th></tr>";
    for (let i = 0; i < memory.length; i += 2) {
        memTable.innerHTML += `<tr><td>${memory[i]}</td><td>${memory[i+1]}</td></tr>`;
    }
}

// Simulation Functions
function runSimulation() {
    // Demo: increase each register by 1
    for (let i = 0; i < NUM_REGISTERS; i++) registers[i] += 1;
    updateTables();
    document.getElementById("status").innerText = "Program executed!";
}

function stepSimulation() {
    // Demo: increase R0 by 1 per step
    registers[0] += 1;
    updateTables();
    document.getElementById("status").innerText = "Step executed!";
}

function resetSimulation() {
    registers = Array(NUM_REGISTERS).fill(0);
    memory = [0,10,1,20,2,30];
    updateTables();
    document.getElementById("status").innerText = "Simulator reset!";
}

// Initial update
updateTables();


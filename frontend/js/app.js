// app.js

let simInitialized = false;

// Utility to create table headers
function createTableHeaders(tableId, headers) {
    const table = document.getElementById(tableId);
    table.innerHTML = "";
    const tr = document.createElement("tr");
    headers.forEach(h => {
        const th = document.createElement("th");
        th.innerText = h;
        tr.appendChild(th);
    });
    table.appendChild(tr);
}

// Utility to fill table data
function fillTable(tableId, data) {
    const table = document.getElementById(tableId);
    // Clear rows except header
    table.querySelectorAll("tr:not(:first-child)").forEach(tr => tr.remove());
    data.forEach(row => {
        const tr = document.createElement("tr");
        row.forEach(cell => {
            const td = document.createElement("td");
            td.innerText = cell;
            tr.appendChild(td);
        });
        table.appendChild(tr);
    });
}

// Initialize simulation with instructions
function initSimulation() {
    const instrText = document.getElementById("instructions").value;
    const instructions = instrText.split("\n").filter(line => line.trim() !== "");

    // Pass instructions to WASM
    const n = instructions.length;
    const ptrs = instructions.map(instr => Module.allocate(Module.intArrayFromString(instr), 'i8', Module.ALLOC_NORMAL));
    const instrArray = Module.allocate(ptrs, 'i32', Module.ALLOC_NORMAL);

    Module.ccall('init_simulation', null, ['number', 'number'], [instrArray, n]);

    simInitialized = true;
    document.getElementById("log").innerText = "Simulation initialized with instructions:\n" + instructions.join("\n");

    // Setup table headers
    createTableHeaders("registerTable", ["Register", "Value"]);
    createTableHeaders("rsTable", ["RS Name", "Busy", "Op", "Vj", "Vk", "Qj", "Qk"]);
    createTableHeaders("robTable", ["ROB Entry", "Busy", "Instruction", "State", "Value"]);
}

// Run next step of simulation
function runNextStep() {
    if (!simInitialized) {
        alert("Please load instructions first!");
        return;
    }

    // Call WASM function to execute next cycle
    Module.ccall('next_cycle', null, [], []);

    // Fetch data from WASM
    const registers = Module.ccall('get_registers', 'array', [], []);
    const rs = Module.ccall('get_reservation_stations', 'array', [], []);
    const rob = Module.ccall('get_rob', 'array', [], []);
    const logText = Module.ccall('get_log', 'string', [], []);

    // Update tables
    fillTable("registerTable", registers);
    fillTable("rsTable", rs);
    fillTable("robTable", rob);
    document.getElementById("log").innerText = logText;
}

// Event listeners
document.getElementById("loadBtn").addEventListener("click", initSimulation);
document.getElementById("runBtn").addEventListener("click", runNextStep);

// When WASM module is ready
Module.onRuntimeInitialized = () => {
    console.log("WASM module loaded and ready!");
};

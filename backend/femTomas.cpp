#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <map>
#include <queue>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <fstream>

using namespace std;

// ==================== Constants and Enums ====================
const int NUM_REGISTERS = 8;
const int ROB_SIZE = 8;
const int MEMORY_SIZE = 65536;

enum InstructionType {
    LOAD, STORE, BEQ, CALL, RET, ADD, SUB, NAND, MUL, INVALID
};

enum ReservationStationState {
    IDLE, BUSY, EXECUTING, WRITING
};

enum ROBState {
    ISSUED, EXECUTING_ROB, WRITTEN, COMMITTED
};

// ==================== Structures ====================
struct Instruction {
    InstructionType type;
    int rA, rB, rC, offset, address;
    string label;
    
    Instruction() : type(INVALID), rA(0), rB(0), rC(0), offset(0), address(0) {}
};
struct InstructionTiming {
    int address;
    int issueCycle;
    int startExecCycle;
    int endExecCycle;
    int writeCycle;
    int commitCycle;
    bool committed;
    
    InstructionTiming() : address(-1), issueCycle(0), startExecCycle(0), 
                          endExecCycle(0), writeCycle(0), commitCycle(0), 
                          committed(false) {}
};


struct RegisterStatus {
    int value;
    int robTag;
    bool valid;
    
    RegisterStatus() : value(0), robTag(-1), valid(true) {}
};

struct ReservationStation {
    int id;
    InstructionType type;
    ReservationStationState state;
    int rA, rB, rC, offset;
    int vj, vk, qj, qk;
    int cyclesRemaining;
    int robTag;
    bool startedExecution;
    
    ReservationStation() : id(-1), type(INVALID), state(IDLE), rA(0), rB(0), rC(0), offset(0),
                          vj(0), vk(0), qj(-1), qk(-1), cyclesRemaining(0), robTag(-1), startedExecution(false) {}
};

struct ROBEntry {
    int tag;
    ROBState state;
    InstructionType type;
      bool isSpeculative;      // Is this instruction speculative?
    int speculativeBranchTag; // Which branch is this speculating on? (-1 if none)
    
    int destination;
    int value;
    bool ready;
    int instructionAddress;
    int issueCycle, startExecCycle, endExecCycle, writeCycle, commitCycle;
    
    ROBEntry() : tag(-1), state(ISSUED), type(INVALID), destination(-1), value(0), 
                ready(false), instructionAddress(0), issueCycle(0), startExecCycle(0),
                endExecCycle(0), writeCycle(0), commitCycle(0),  isSpeculative(false), speculativeBranchTag(-1) {}
};

// ==================== Tomasulo Simulator ====================
class TomasuloSimulator {
private:
    vector<Instruction> instructionQueue;
    vector<ReservationStation> reservationStations;
    vector<ROBEntry> reorderBuffer;
    vector<RegisterStatus> registerFile;
    vector<int> memory;
    int lastUnresolvedBranchROB; 
        int nextInstructionToIssue; 
  bool branchMispredicted;
    int branchTarget;
 vector<InstructionTiming> instructionTimings; 
    int currentCycle;
    int instructionsCompleted;
    int conditionalBranches;
    int branchMispredictions;
    int totalCycles;
    
    int programCounter;
    bool stall;
    bool debugMode;
    
    map<string, InstructionType> instructionMap;
    
public:
    TomasuloSimulator() {
        initialize();
    }
    
    void initialize() {
        memory.resize(MEMORY_SIZE, 0);
        lastUnresolvedBranchROB = -1;

        registerFile.resize(NUM_REGISTERS);
        for (int i = 0; i < NUM_REGISTERS; i++) {
            registerFile[i].value = 0;
            registerFile[i].valid = true;
            registerFile[i].robTag = -1;
        }
        registerFile[0].value = 0;  // R0 always 0
        
        initializeInstructionMap();
        initializeReservationStations();
        
        reorderBuffer.resize(ROB_SIZE);
        for (int i = 0; i < ROB_SIZE; i++) {
            reorderBuffer[i].tag = i;
            reorderBuffer[i].state = ISSUED;
        }
        
        currentCycle = 0;
        instructionsCompleted = 0;
        conditionalBranches = 0;
        branchMispredictions = 0;
        totalCycles = 0;
        
        programCounter = 0;
        stall = false;
        debugMode = true;
    }
    
    void initializeInstructionMap() {
        instructionMap["LOAD"] = LOAD;
        instructionMap["STORE"] = STORE;
        instructionMap["BEQ"] = BEQ;
        instructionMap["CALL"] = CALL;
        instructionMap["RET"] = RET;
        instructionMap["ADD"] = ADD;
        instructionMap["SUB"] = SUB;
        instructionMap["NAND"] = NAND;
        instructionMap["MUL"] = MUL;
    }
    
    void initializeReservationStations() {
        reservationStations.clear();
        
        // 2 LOAD units
        for (int i = 0; i < 2; i++) {
            ReservationStation rs;
            rs.id = reservationStations.size();
            rs.type = LOAD;
            reservationStations.push_back(rs);
        }
        
        // 1 STORE unit
        ReservationStation storeRS;
        storeRS.id = reservationStations.size();
        storeRS.type = STORE;
        reservationStations.push_back(storeRS);
        
        // 2 BEQ units
        for (int i = 0; i < 2; i++) {
            ReservationStation rs;
            rs.id = reservationStations.size();
            rs.type = BEQ;
            reservationStations.push_back(rs);
        }
        
        // 1 CALL/RET unit
        ReservationStation callRetRS;
        callRetRS.id = reservationStations.size();
        callRetRS.type = CALL;
        reservationStations.push_back(callRetRS);
        
        // 4 ADD/SUB units
        for (int i = 0; i < 4; i++) {
            ReservationStation rs;
            rs.id = reservationStations.size();
            rs.type = ADD;
            reservationStations.push_back(rs);
        }
        
        // 2 NAND units
        for (int i = 0; i < 2; i++) {
            ReservationStation rs;
            rs.id = reservationStations.size();
            rs.type = NAND;
            reservationStations.push_back(rs);
        }
        
        // 1 MUL unit
        ReservationStation mulRS;
        mulRS.id = reservationStations.size();
        mulRS.type = MUL;
        reservationStations.push_back(mulRS);
    }
    
    Instruction parseInstruction(const string& line, int address) {
        Instruction instr;
        instr.address = address;
        
        stringstream ss(line);
        string opcode;
        ss >> opcode;
        
        for (char &c : opcode) c = toupper(c);
        
        if (instructionMap.find(opcode) != instructionMap.end()) {
            instr.type = instructionMap[opcode];
        } else {
            instr.type = INVALID;
            return instr;
        }
        
        switch (instr.type) {
            case LOAD: {
                string reg, offsetReg;
                ss >> reg;
                if (!reg.empty() && reg.back() == ',') reg.pop_back();
                instr.rA = reg[1] - '0';
                
                ss >> offsetReg;
                size_t paren = offsetReg.find('(');
                if (paren != string::npos) {
                    instr.offset = stoi(offsetReg.substr(0, paren));
                    size_t rstart = offsetReg.find('R', paren);
                    if (rstart != string::npos) {
                        instr.rB = offsetReg[rstart + 1] - '0';
                    }
                }
                break;
            }
            case STORE: {
                string reg, offsetReg;
                ss >> reg;
                if (!reg.empty() && reg.back() == ',') reg.pop_back();
                instr.rA = reg[1] - '0';
                
                ss >> offsetReg;
                size_t paren = offsetReg.find('(');
                if (paren != string::npos) {
                    instr.offset = stoi(offsetReg.substr(0, paren));
                    size_t rstart = offsetReg.find('R', paren);
                    if (rstart != string::npos) {
                        instr.rB = offsetReg[rstart + 1] - '0';
                    }
                }
                break;
            }
            case BEQ: {
                string reg1, reg2;
                ss >> reg1 >> reg2;
                if (!reg1.empty() && reg1.back() == ',') reg1.pop_back();
                if (!reg2.empty() && reg2.back() == ',') reg2.pop_back();
                instr.rA = reg1[1] - '0';
                instr.rB = reg2[1] - '0';
                
                int offset;
                ss >> offset;
                instr.offset = offset;
                break;
            }
            case CALL: {
                int label;
                ss >> label;
                instr.offset = label;
                break;
            }
            case RET:
                break;
            case ADD:
            case SUB:
            case NAND:
            case MUL: {
                string reg1, reg2, reg3;
                ss >> reg1 >> reg2 >> reg3;
                if (!reg1.empty() && reg1.back() == ',') reg1.pop_back();
                if (!reg2.empty() && reg2.back() == ',') reg2.pop_back();
                instr.rA = reg1[1] - '0';
                instr.rB = reg2[1] - '0';
                instr.rC = reg3[1] - '0';
                break;
            }
            default:
                break;
        }
        
        stringstream labelSS;
        labelSS << opcode << " ";
        switch (instr.type) {
            case LOAD:
                labelSS << "R" << instr.rA << ", " << instr.offset << "(R" << instr.rB << ")";
                break;
            case STORE:
                labelSS << "R" << instr.rA << ", " << instr.offset << "(R" << instr.rB << ")";
                break;
            case BEQ:
                labelSS << "R" << instr.rA << ", R" << instr.rB << ", " << instr.offset;
                break;
            case CALL:
                labelSS << instr.offset;
                break;
            case RET:
                labelSS << "";
                break;
            case ADD:
            case SUB:
            case NAND:
            case MUL:
                labelSS << "R" << instr.rA << ", R" << instr.rB << ", R" << instr.rC;
                break;
            default:
                break;
        }
        instr.label = labelSS.str();
        
        return instr;
    }
    
    int allocateROBEntry() {
        for (int i = 0; i < ROB_SIZE; i++) {
            if (reorderBuffer[i].state == ISSUED) {
                reorderBuffer[i].state = EXECUTING_ROB;
                return i;
            }
        }
        return -1;
    }
    

int findFreeReservationStation(InstructionType type) {
    for (auto& rs : reservationStations) {
        if (rs.state == IDLE) {
            bool canHandle = false;
            
            switch (type) {
                case LOAD:
                    canHandle = (rs.type == LOAD);
                    break;
                case STORE:
                    canHandle = (rs.type == STORE);
                    break;
                case BEQ:
                    canHandle = (rs.type == BEQ);
                    break;
                case CALL:
                case RET:
                    canHandle = (rs.type == CALL);
                    break;
                case ADD:
                case SUB:
                    canHandle = (rs.type == ADD);
                    break;
                case NAND:
                    canHandle = (rs.type == NAND);
                    break;
                case MUL:
                    canHandle = (rs.type == MUL);
                    break;
                default:
                    break;
            }
            
            if (canHandle) {
                return rs.id;
            }
        }
    }
    return -1;
}


    bool issueInstruction(Instruction& instr) {
        int robTag = allocateROBEntry();
        if (robTag == -1) {
            if (debugMode) cout << "  ROB full\n";
            return false;
        }
        
        int rsId = findFreeReservationStation(instr.type);
        if (rsId == -1) {
            reorderBuffer[robTag].state = ISSUED;
            if (debugMode) cout << "  No RS available\n";
            return false;
        }
        
        ReservationStation& rs = reservationStations[rsId];
        ROBEntry& robEntry = reorderBuffer[robTag];
        
        if (debugMode) {
            cout << "Cycle " << currentCycle << ": Issue " << instr.label 
                 << " -> RS" << rsId << ", ROB" << robTag << "\n";
        }
        
        // Configure RS
        rs.state = BUSY;
        rs.type = instr.type; 
        rs.rA = instr.rA;
        rs.rB = instr.rB;
        rs.rC = instr.rC;
        rs.offset = instr.offset;
        rs.robTag = robTag;
        rs.cyclesRemaining = getCyclesForType(instr.type);
        rs.startedExecution = false;
        
        // Configure ROB
        robEntry.type = instr.type;
        robEntry.instructionAddress = instr.address;
        robEntry.issueCycle = currentCycle;
         instructionTimings[instr.address].issueCycle = currentCycle;
        robEntry.destination = -1;
        

         // NEW: Track speculation
    if (instr.type == BEQ) {
        // This is a branch - mark it as the last unresolved branch
        lastUnresolvedBranchROB = robTag;
        robEntry.isSpeculative = false;
        robEntry.speculativeBranchTag = -1;
    } else if (lastUnresolvedBranchROB != -1) {
        // There's an unresolved branch, so this instruction is speculative
        robEntry.isSpeculative = true;
        robEntry.speculativeBranchTag = lastUnresolvedBranchROB;
    } else {
        robEntry.isSpeculative = false;
        robEntry.speculativeBranchTag = -1;
    }
        // Handle dependencies
        if (instr.type == LOAD) {
            // LOAD rA, offset(rB)
            if (registerFile[instr.rB].valid) {
                rs.vj = registerFile[instr.rB].value;
                rs.qj = -1;
            } else {
                rs.qj = registerFile[instr.rB].robTag;
            }
            robEntry.destination = instr.rA;
            if (instr.rA > 0) {
                registerFile[instr.rA].robTag = robTag;
                registerFile[instr.rA].valid = false;
            }
        } 
        else if (instr.type == STORE) {
            // STORE rA, offset(rB)
            if (registerFile[instr.rB].valid) {
                rs.vj = registerFile[instr.rB].value;
                rs.qj = -1;
            } else {
                rs.qj = registerFile[instr.rB].robTag;
            }
            if (registerFile[instr.rA].valid) {
                rs.vk = registerFile[instr.rA].value;
                rs.qk = -1;
            } else {
                rs.qk = registerFile[instr.rA].robTag;
            }
        }
        else if (instr.type >= ADD && instr.type <= MUL) {
            // ALU: rA = rB op rC
            if (registerFile[instr.rB].valid) {
                rs.vj = registerFile[instr.rB].value;
                rs.qj = -1;
            } else {
                rs.qj = registerFile[instr.rB].robTag;
            }
            if (registerFile[instr.rC].valid) {
                rs.vk = registerFile[instr.rC].value;
                rs.qk = -1;
            } else {
                rs.qk = registerFile[instr.rC].robTag;
            }
            robEntry.destination = instr.rA;
            if (instr.rA > 0) {
                registerFile[instr.rA].robTag = robTag;
                registerFile[instr.rA].valid = false;
            }
        }
        else if (instr.type == BEQ) {
            if (registerFile[instr.rA].valid) {
                rs.vj = registerFile[instr.rA].value;
                rs.qj = -1;
            } else {
                rs.qj = registerFile[instr.rA].robTag;
            }
            if (registerFile[instr.rB].valid) {
                rs.vk = registerFile[instr.rB].value;
                rs.qk = -1;
            } else {
                rs.qk = registerFile[instr.rB].robTag;
            }
            
            conditionalBranches++;
        }
        
        return true;
    }

    void flushSpeculativeInstructions(int branchROBTag) {
    if (debugMode) {
        cout << "  FLUSHING speculative instructions after branch ROB" << branchROBTag << "\n";
    }
    
    // Flush all ROB entries that are speculative on this branch
    for (int i = 0; i < ROB_SIZE; i++) {
        ROBEntry& entry = reorderBuffer[i];
        if (entry.speculativeBranchTag == branchROBTag && entry.type != INVALID) {
            if (debugMode) {
                cout << "    Flush ROB" << i << " (" << instructionTypeToString(entry.type) << ")\n";
            }
            
            // Free reserved register
            // Since registers are only updated at commit (not broadcast),
            // we just need to mark the register as available again
            if (entry.destination >= 0 && entry.destination < NUM_REGISTERS) {
                if (registerFile[entry.destination].robTag == i) {
                    registerFile[entry.destination].valid = true;
                    registerFile[entry.destination].robTag = -1;
                }
            }
            
            // Clear the ROB entry completely
            entry.state = ISSUED;
            entry.type = INVALID;
            entry.ready = false;
            entry.isSpeculative = false;
            entry.speculativeBranchTag = -1;
            entry.destination = -1;
        }
    }
    
    // Flush all RS entries for speculative instructions
    for (auto& rs : reservationStations) {
        if (rs.robTag != -1) {
            ROBEntry& entry = reorderBuffer[rs.robTag];
            if (entry.speculativeBranchTag == branchROBTag || entry.type == INVALID) {
                if (debugMode) {
                    cout << "    Flush RS" << rs.id << "\n";
                }
                rs.state = IDLE;
                rs.robTag = -1;
                rs.startedExecution = false;
                rs.qj = -1;
                rs.qk = -1;
            }
        }
    }
}

void executeStage() {
    for (auto& rs : reservationStations) {
        // 1. Check if READY to start execution
        if (rs.state == BUSY && rs.qj == -1 && rs.qk == -1 && !rs.startedExecution) {
            rs.state = EXECUTING;
            rs.startedExecution = true;
            reorderBuffer[rs.robTag].startExecCycle = currentCycle;
            instructionTimings[reorderBuffer[rs.robTag].instructionAddress].startExecCycle = currentCycle;
            if (debugMode) {
                cout << "Cycle " << currentCycle << ": Start exec RS" << rs.id 
                     << " (ROB" << rs.robTag << "), cycles=" << rs.cyclesRemaining << "\n";
            }
            continue;
        }
        
        // 2. Continue execution for already executing RS
        if (rs.state == EXECUTING && rs.cyclesRemaining > 0) {
            rs.cyclesRemaining--;
            
            if (debugMode && rs.cyclesRemaining > 0) {
                cout << "  RS" << rs.id << " exec, " << rs.cyclesRemaining << " cycles left\n";
            }
            
            if (rs.cyclesRemaining == 0) {
                int result = executeInstruction(rs);
                ROBEntry& entry = reorderBuffer[rs.robTag];
                entry.value = result;
                entry.endExecCycle = currentCycle;
                   instructionTimings[entry.instructionAddress].endExecCycle = currentCycle;
                entry.ready = true;
                rs.state = WRITING;
                
                if (debugMode) {
                    cout << "Cycle " << currentCycle << ": Finish exec RS" << rs.id 
                         << ", result=" << result << "\n";
                }
            }
        }
        else if (rs.state == BUSY && debugMode) {
            cout << "  RS" << rs.id << " waiting: qj=" << rs.qj << ", qk=" << rs.qk << "\n";
        }
    }
}

    int executeInstruction(ReservationStation& rs) {
        int result = 0;
        
        switch (rs.type) {
            case LOAD: {
                int address = rs.vj + rs.offset;
                result = memory[address];
                if (debugMode) {
                    cout << "    LOAD: " << rs.vj << " + " << rs.offset << " = addr " 
                         << address << " -> " << result << "\n";
                }
                break;
            }
            case STORE: {
                int address = rs.vj + rs.offset;
                result = rs.vk;
                reorderBuffer[rs.robTag].destination = address;
                if (debugMode) {
                    cout << "    STORE: addr=" << address << ", value=" << result << "\n";
                }
                break;
            }
            case BEQ: {
                result = (rs.vj == rs.vk) ? 1 : 0;
                if (result == 1) {
                    branchMispredictions++;
                    if (debugMode) cout << "    BEQ misprediction!\n";
                }
                break;
            }
            case ADD:
                result = rs.vj + rs.vk;
                if (debugMode) cout << "    ADD: " << rs.vj << " + " << rs.vk << " = " << result << "\n";
                break;
            case SUB:
                result = rs.vj - rs.vk;
                break;
            case NAND:
                result = ~(rs.vj & rs.vk);
                break;
            case MUL: {
                int32_t product = (int32_t)rs.vj * (int32_t)rs.vk;
                result = product & 0xFFFF;
                break;
            }
            default:
                break;
        }
        
        return result & 0xFFFF;
    }
    
    int getCyclesForType(InstructionType type) {
        switch (type) {
            case LOAD: return 6;
            case STORE: return 6;
            case BEQ: return 1;
            case CALL:
            case RET: return 1;
            case ADD:
            case SUB: return 2;
            case NAND: return 1;
            case MUL: return 12;
            default: return 1;
        }
    }
    

    void writeResultStage() {
    for (auto& rs : reservationStations) {
        if (rs.state == WRITING) {
            ROBEntry& robEntry = reorderBuffer[rs.robTag];
            
            // Check if this instruction should be flushed
            // (its branch predecessor mispredicted)
            if (robEntry.isSpeculative && robEntry.speculativeBranchTag != -1) {
                // Check if the branch this depends on has mispredicted
                ROBEntry& branchEntry = reorderBuffer[robEntry.speculativeBranchTag];
                if (branchEntry.type == BEQ && branchEntry.ready && branchEntry.value == 1) {
                    // Branch mispredicted - don't broadcast, just mark for flush
                    if (debugMode) {
                        cout << "Cycle " << currentCycle << ": SKIP write RS" << rs.id 
                             << " (ROB" << rs.robTag << ") - speculative on mispredicted branch\n";
                    }
                    // Clear RS but don't broadcast
                    rs.state = IDLE;
                    rs.robTag = -1;
                    rs.startedExecution = false;
                    continue;
                }
            }
            
            // Normal broadcast
            broadcastResult(rs.robTag, robEntry.value);
            robEntry.writeCycle = currentCycle;
             instructionTimings[robEntry.instructionAddress].writeCycle = currentCycle; 
            if (debugMode) {
                cout << "Cycle " << currentCycle << ": Write result RS" << rs.id 
                     << " (ROB" << rs.robTag << ") = " << robEntry.value << "\n";
            }
            
            // Clear RS
            rs.state = IDLE;
            rs.robTag = -1;
            rs.startedExecution = false;
        }
    }
}

void broadcastResult(int robTag, int value) {
    // Update reservation stations ONLY
    for (auto& rs : reservationStations) {
        if (rs.qj == robTag) {
            rs.vj = value;
            rs.qj = -1;
            if (debugMode) {
                cout << "    Broadcast to RS" << rs.id << ": vj = " << value << "\n";
            }
        }
        if (rs.qk == robTag) {
            rs.vk = value;
            rs.qk = -1;
            if (debugMode) {
                cout << "    Broadcast to RS" << rs.id << ": vk = " << value << "\n";
            }
        }
    }
}

void simulate() {
    nextInstructionToIssue = 0;
    branchMispredicted = false;
    branchTarget = -1;
    int maxInstructions = instructionQueue.size();
    
    cout << "\n=== Starting Simulation ===\n";
    cout << "Total instructions: " << maxInstructions << "\n";

    bool allWorkDone = false;
    
    while (!allWorkDone) {
        currentCycle++;

        if (debugMode) {
            cout << "\n--- Cycle " << currentCycle << " ---\n";
        }

        commitStage();  // May set branchMispredicted=true and branchTarget
        writeResultStage();
        executeStage();

        // Handle branch misprediction redirect
        if (branchMispredicted) {
            nextInstructionToIssue = branchTarget;
            branchMispredicted = false;
            if (debugMode) {
                cout << "  Redirecting issue to address " << branchTarget << "\n";
            }
        }

        // Issue next instruction
        if (nextInstructionToIssue < maxInstructions) {
            if (issueInstruction(instructionQueue[nextInstructionToIssue])) {
                nextInstructionToIssue++;
            }
        }

        // Check termination
        if (nextInstructionToIssue >= maxInstructions) {
            bool allROBFree = true;
            for (const auto& entry : reorderBuffer) {
                if (entry.type != INVALID && entry.state != ISSUED) {
                    allROBFree = false;
                    break;
                }
            }
            
            bool allRSIdle = true;
            for (const auto& rs : reservationStations) {
                if (rs.state != IDLE) {
                    allRSIdle = false;
                    break;
                }
            }
            
            if (allROBFree && allRSIdle) {
                allWorkDone = true;
            }
        }

        totalCycles = currentCycle;

        if (currentCycle > 100) {  // Increase safety limit
            cout << "Warning: Stopping after 100 cycles\n";
            break;
        }
    }
    
    cout << "\n=== Simulation Complete ===\n";
    cout << "Total cycles: " << totalCycles << "\n";
}

    
void commitStage() {
    for (int i = 0; i < ROB_SIZE; i++) {
        ROBEntry& entry = reorderBuffer[i];
        if (entry.type != INVALID && entry.state == EXECUTING_ROB && 
            entry.ready && entry.writeCycle < currentCycle) {
            
            // Check for branch misprediction FIRST
            if (entry.type == BEQ && entry.value == 1) {
                if (debugMode) {
                    cout << "Cycle " << currentCycle << ": Commit ROB" << i
                         << " (BEQ MISPREDICTION - FLUSHING)\n";
                }
                
                // Flush all speculative instructions
                flushSpeculativeInstructions(i);
                
                // Calculate branch target and set redirection
                int branchPC = entry.instructionAddress;
                int branchOffset = 0;
                
                // Find the original BEQ instruction to get its offset
                for (const auto& instr : instructionQueue) {
                    if (instr.address == branchPC && instr.type == BEQ) {
                        branchOffset = instr.offset;
                        break;
                    }
                }
                
                // BEQ target: PC + 1 + offset
                int targetAddress = branchPC + 1 + branchOffset;
                
                if (debugMode) {
                    cout << "  Branch taken: PC=" << branchPC 
                         << " + 1 + offset(" << branchOffset 
                         << ") = target address " << targetAddress << "\n";
                }
                
                // Set redirection flags for next cycle
                branchMispredicted = true;
                branchTarget = targetAddress;
                
                // Commit the branch itself - update timing
                entry.commitCycle = currentCycle;
                instructionTimings[entry.instructionAddress].commitCycle = currentCycle;
                instructionTimings[entry.instructionAddress].committed = true;
                instructionsCompleted++;
                
                // Clear the branch from ROB
                if (lastUnresolvedBranchROB == i) {
                    lastUnresolvedBranchROB = -1;
                }
                
                // Free ROB entry
                entry.state = ISSUED;
                entry.type = INVALID;
                entry.ready = false;
                
                break;  // One commit per cycle
            }
            
            // SKIP committing if this instruction is speculative on a mispredicted branch
            if (entry.isSpeculative && entry.speculativeBranchTag != -1) {
                ROBEntry& branchEntry = reorderBuffer[entry.speculativeBranchTag];
                if (branchEntry.type == BEQ && branchEntry.ready && branchEntry.value == 1) {
                    // This instruction is speculative on a mispredicted branch
                    // It should have been flushed - skip it
                    if (debugMode) {
                        cout << "Cycle " << currentCycle << ": SKIP commit ROB" << i
                             << " (speculative on mispredicted branch)\n";
                    }
                    // Don't increment instructionsCompleted
                    // Just clear it
                    entry.state = ISSUED;
                    entry.type = INVALID;
                    entry.ready = false;
                    entry.isSpeculative = false;
                    entry.speculativeBranchTag = -1;
                    continue; // Try next ROB entry
                }
            }
            
            // Normal commit
            if (debugMode) {
                cout << "Cycle " << currentCycle << ": Commit ROB" << i
                     << " (" << instructionTypeToString(entry.type) << ")\n";
            }
            
            // **UPDATE ARCHITECTURAL STATE AT COMMIT TIME**
            // Update register file 
            if (entry.destination >= 0 && entry.destination < NUM_REGISTERS) {
                if (registerFile[entry.destination].robTag == i) {
                    registerFile[entry.destination].value = entry.value;
                    registerFile[entry.destination].valid = true;
                    registerFile[entry.destination].robTag = -1;
                    if (debugMode) {
                        cout << "    Update R" << entry.destination 
                             << " = " << entry.value << "\n";
                    }
                }
            }
            
            // Handle STORE - write to memory at commit time
            if (entry.type == STORE) {
                if (entry.destination >= 0 && entry.destination < MEMORY_SIZE) {
                    memory[entry.destination] = entry.value;
                    if (debugMode) {
                        cout << "    Store to memory[" << entry.destination 
                             << "] = " << entry.value << "\n";
                    }
                }
            }
            
            // Mark committed and update timing
            entry.commitCycle = currentCycle;
            instructionTimings[entry.instructionAddress].commitCycle = currentCycle;
            instructionTimings[entry.instructionAddress].committed = true;
            instructionsCompleted++;
            
            // Clear unresolved branch tracking if this was a branch
            if (entry.type == BEQ && lastUnresolvedBranchROB == i) {
                lastUnresolvedBranchROB = -1;
            }
            
            // Free ROB entry
            entry.state = ISSUED;
            entry.type = INVALID;
            entry.ready = false;
            entry.destination = -1;
            
            break;  // One commit per cycle
        }
    }
}
    
    string instructionTypeToString(InstructionType type) {
        switch (type) {
            case LOAD: return "LOAD";
            case STORE: return "STORE";
            case BEQ: return "BEQ";
            case CALL: return "CALL";
            case RET: return "RET";
            case ADD: return "ADD";
            case SUB: return "SUB";
            case NAND: return "NAND";
            case MUL: return "MUL";
            default: return "INVALID";
        }
    }
    
    
    void loadDefaultProgram() {
    cout << "Loading default test program...\n";
    
    instructionQueue.clear();
    instructionQueue.push_back(parseInstruction("LOAD R1, 0(R0)", 0));
    instructionQueue.push_back(parseInstruction("ADD R2, R1, R1", 1));
    instructionQueue.push_back(parseInstruction("STORE R2, 4(R0)", 2));
    
    instructionTimings.clear();
    instructionTimings.resize(instructionQueue.size());
    for (int i = 0; i < (int)instructionQueue.size(); i++) {
        instructionTimings[i].address = i;
    }
    
    // Initialize memory
    memory[0] = 10;
    
    // Initialize registers
    for (int i = 0; i < NUM_REGISTERS; i++) {
        registerFile[i].value = 0;
        registerFile[i].valid = true;
        registerFile[i].robTag = -1;
    }
    registerFile[0].value = 0;  // R0 always 0
    
    cout << "Default program loaded:\n";
    for (int i = 0; i < instructionQueue.size(); i++) {
        cout << i << ": " << instructionQueue[i].label << "\n";
    }
    
    cout << "Memory[0] = 10\n";
}


void loadProgramFromFile() {
    cout << "Loading program from program.txt and memory.txt...\n";

    instructionQueue.clear();

    // 1) Load instructions
    ifstream progFile("program.txt");
    if (!progFile.is_open()) {
        cout << "Error: Cannot open program.txt\n";
        return;
    }

    string line;
    int address = 0;
    while (getline(progFile, line)) {
        // ignore empty / whitespace-only lines
        bool empty = true;
        for (char c : line) {
            if (!isspace(static_cast<unsigned char>(c))) { empty = false; break; }
        }
        if (!empty) {
            instructionQueue.push_back(parseInstruction(line, address++));
        }
    }
    progFile.close();

    instructionTimings.clear();
    instructionTimings.resize(instructionQueue.size());
    for (int i = 0; i < (int)instructionQueue.size(); i++) {
        instructionTimings[i].address = i;
    }

    // 2) Clear and load memory initial contents
    for (int i = 0; i < MEMORY_SIZE; i++) {
        memory[i] = 0;
    }

    ifstream memFile("memory.txt");
    if (!memFile.is_open()) {
        cout << "Warning: Cannot open memory.txt (starting with all zeros)\n";
    } else {
        int addr, value;
        while (memFile >> addr >> value) {
            if (addr >= 0 && addr < MEMORY_SIZE) {
                memory[addr] = value;
            }
        }
        memFile.close();
    }

    // 3) Reset registers
    for (int i = 0; i < NUM_REGISTERS; i++) {
        registerFile[i].value = 0;
        registerFile[i].valid = true;
        registerFile[i].robTag = -1;
    }
    registerFile[0].value = 0; // R0 always 0

    cout << "Program loaded from files:\n";
    for (int i = 0; i < (int)instructionQueue.size(); i++) {
        cout << i << ": " << instructionQueue[i].label << "\n";
    }
}

    
    void printResults() {
    cout << "\n=== SIMULATION RESULTS ===\n\n";
    
    cout << "Instruction Timeline:\n";
    cout << "================================================================================\n";
    cout << left << setw(5) << "Addr" << setw(20) << "Instruction" 
         << setw(8) << "Issue" << setw(8) << "Exec" 
         << setw(8) << "Finish" << setw(8) << "Write" 
         << setw(8) << "Commit" << "\n";
    cout << "================================================================================\n";
    
    for (const auto& instr : instructionQueue) {
        const InstructionTiming& timing = instructionTimings[instr.address];
        
        cout << left << setw(5) << instr.address 
             << setw(20) << instr.label 
             << setw(8) << (timing.issueCycle > 0 ? to_string(timing.issueCycle) : "-")
             << setw(8) << (timing.startExecCycle > 0 ? to_string(timing.startExecCycle) : "-")
             << setw(8) << (timing.endExecCycle > 0 ? to_string(timing.endExecCycle) : "-")
             << setw(8) << (timing.writeCycle > 0 ? to_string(timing.writeCycle) : "-")
             << setw(8) << (timing.commitCycle > 0 ? to_string(timing.commitCycle) : "-") << "\n";
    }
    
    // ... rest of printResults() stays the same ...
    
    cout << "\n=== Performance Metrics ===\n";
    cout << "Total Cycles: " << totalCycles << "\n";
    cout << "Completed instructions: " << instructionsCompleted << "\n";
    
    double ipc = (instructionsCompleted > 0 && totalCycles > 0) ? 
                 (double)instructionsCompleted / totalCycles : 0.0;
    cout << "IPC: " << fixed << setprecision(3) << ipc << "\n";
    
    cout << "Conditional Branches: " << conditionalBranches << "\n";
    cout << "Branch Mispredictions: " << branchMispredictions << "\n";
    
    double mispredictionRate = (conditionalBranches > 0) ?
                              (double)branchMispredictions / conditionalBranches * 100 : 0.0;
    cout << "Branch Misprediction Rate: " << fixed << setprecision(2) 
         << mispredictionRate << "%\n";
    
    cout << "\nRegisters:\n";
    for (int i = 0; i < NUM_REGISTERS; i++) {
        cout << "R" << i << ": " << registerFile[i].value;
        if (!registerFile[i].valid) {
            cout << " (Waiting for ROB " << registerFile[i].robTag << ")";
        }
        cout << "\n";
    }
    
    cout << "\nMemory[0-9]:\n";
    for (int i = 0; i < 10; i++) {
        cout << "M[" << i << "]: " << memory[i] << "\n";
    }
}

    void run() {
        char choice;
        cout << "Load default test program? (y/n): ";
        cin >> choice;
        cin.ignore();
        
           if (choice == 'y' || choice == 'Y') {
        loadDefaultProgram();
    } else {
        loadProgramFromFile();
    }

        simulate();
        printResults();
    }
};

// ==================== Main Function ====================
int main() {
    cout << "CSCE 3301 - Computer Architecture - Fall 2025\n";
    cout << "Project 2: femTomas - Tomasulo Algorithm Simulator\n";
    cout << "==============================================\n\n";
    
    TomasuloSimulator simulator;
    simulator.run();
 
    cout << "\n\nExecution finished. Press Enter to exit...";
cin.ignore();
cin.get();

    return 0;
}
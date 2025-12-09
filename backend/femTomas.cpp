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
    bool isSpeculative;
    int speculativeBranchTag;
    
    int destination;
    int value;
    bool ready;
    int instructionAddress;
    int issueCycle, startExecCycle, endExecCycle, writeCycle, commitCycle;
    
    // For CALL instructions, we need to store the offset
    int callOffset;
    
    ROBEntry() : tag(-1), state(ISSUED), type(INVALID), destination(-1), value(0), 
                ready(false), instructionAddress(0), issueCycle(0), startExecCycle(0),
                endExecCycle(0), writeCycle(0), commitCycle(0), isSpeculative(false), 
                speculativeBranchTag(-1), callOffset(0) {}
};

// ==================== Tomasulo Simulator ====================
class TomasuloSimulator {
private:
int robHead;
int robTail;
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
    int programStartAddress;
    int commitCyclesRemaining = 0;
int committingROBIndex = -1;

    int programCounter;
    bool stall;
    bool debugMode;
    
    map<string, InstructionType> instructionMap;
    bool isROBEmpty() const {
    for (const auto& entry : reorderBuffer) {
        if (entry.type != INVALID) {
            return false;
        }
    }
    return true;
}
    // Helper function to sign-extend 5-bit offset
    int signExtend5To16(int offset) {
        offset &= 0x1F;  // Keep only lower 5 bits
        if (offset & 0x10) {  // If negative (bit 4 = 1)
            return offset | 0xFFE0;  // Sign extend to 16 bits
        }
        return offset;
    }
    
    // Helper function to sign-extend 7-bit offset
    int signExtend7To16(int offset) {
        offset &= 0x7F;  // Keep only lower 7 bits
        if (offset & 0x40) {  // If negative (bit 6 = 1)
            return offset | 0xFF80;  // Sign extend to 16 bits
        }
        return offset;
    }
    
    // Helper function for 16-bit arithmetic
    int to16Bit(int value) {
        return value & 0xFFFF;
    }
    
public:
    TomasuloSimulator() {
        initialize();
    }
    
   void initialize() {
    memory.resize(MEMORY_SIZE, 0);
    lastUnresolvedBranchROB = -1;
    programStartAddress = 0;

    registerFile.resize(NUM_REGISTERS);
    for (int i = 0; i < NUM_REGISTERS; i++) {
        registerFile[i].value = 0;
        registerFile[i].valid = true;
        registerFile[i].robTag = -1;
    }
    registerFile[0].value = 0;  // R0 always 0
    
    robHead = 0;
    robTail = 0;
    
    initializeInstructionMap();
    initializeReservationStations();
    
    reorderBuffer.resize(ROB_SIZE);
    for (int i = 0; i < ROB_SIZE; i++) {
        reorderBuffer[i].tag = i;
        reorderBuffer[i].state = ISSUED;
        reorderBuffer[i].type = INVALID;    // CRITICAL: Start as INVALID
        reorderBuffer[i].ready = false;     // CRITICAL: Start as not ready
        reorderBuffer[i].isSpeculative = false;
        reorderBuffer[i].speculativeBranchTag = -1;
    }
    
    currentCycle = 0;
    instructionsCompleted = 0;
    conditionalBranches = 0;
    branchMispredictions = 0;
    totalCycles = 0;
    
    programCounter = 0;
    stall = false;
    debugMode = true;
    
    commitCyclesRemaining = 0;
    committingROBIndex = -1;
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
                // instr.rA = reg[1] - '0';
                instr.rA = stoi(reg.substr(1));

                ss >> offsetReg;
                size_t paren = offsetReg.find('(');
                if (paren != string::npos) {
                    string offsetStr = offsetReg.substr(0, paren);
                    instr.offset = stoi(offsetStr);
                    instr.offset = signExtend5To16(instr.offset);  // 5-bit signed
                    
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
                    string offsetStr = offsetReg.substr(0, paren);
                    instr.offset = stoi(offsetStr);
                    instr.offset = signExtend5To16(instr.offset);  // 5-bit signed
                    
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
                if (ss >> offset) {
                    instr.offset = signExtend5To16(offset);  // 5-bit signed
                }
                break;
            }
            case CALL: {
                int label;
                if (ss >> label) {
                    instr.offset = signExtend7To16(label);  // 7-bit signed
                }
                break;
            }
            case RET:
                // No operands
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
    // Check if ROB is full
    if ((robTail + 1) % ROB_SIZE == robHead && 
        reorderBuffer[robHead].type != INVALID) {
        return -1;  // ROB full
    }
    
    // Allocate at tail
    int allocated = robTail;
    robTail = (robTail + 1) % ROB_SIZE;
    return allocated;
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
        
        // Store CALL offset in ROB
        if (instr.type == CALL) {
            robEntry.callOffset = instr.offset;
        }
        
        // Track speculation
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
        if (instr.type == BEQ) {
    lastUnresolvedBranchROB = robTag;
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
        else if (instr.type == CALL) {
            // CALL label - uses R1 for return address
            robEntry.destination = 1;  // R1
            registerFile[1].robTag = robTag;
            registerFile[1].valid = false;
        }
        else if (instr.type == RET) {
            // RET - depends on R1
            if (registerFile[1].valid) {
                rs.vj = registerFile[1].value;
                rs.qj = -1;
            } else {
                rs.qj = registerFile[1].robTag;
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
            if (entry.destination >= 0 && entry.destination < NUM_REGISTERS) {
                if (registerFile[entry.destination].robTag == i) {
                    registerFile[entry.destination].valid = true;
                    registerFile[entry.destination].robTag = -1;
                }
            }
            
            // Clear the ROB entry
            entry.state = ISSUED;
            entry.type = INVALID;
            entry.ready = false;
            entry.isSpeculative = false;
            entry.speculativeBranchTag = -1;
            entry.destination = -1;
        }
    }
    
    // Reset tail to one past the branch (circular)
    robTail = (branchROBTag + 1) % ROB_SIZE;  // Add this line
    
    // Flush matching RS entries
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
                cout << "Cycle " << currentCycle << ": Start exec RS " << rs.id 
                     << " ROB " << rs.robTag << ", cycles=" << rs.cyclesRemaining << "\n";
            }
        }

        // 2. Continue execution for multi-cycle instructions
        if (rs.state == EXECUTING && rs.cyclesRemaining > 0) {
            if (debugMode && rs.cyclesRemaining > 0) {
                cout << "RS " << rs.id << " exec, " << rs.cyclesRemaining << " cycles left\n";
            }
            
            rs.cyclesRemaining--;
            
            // Check if execution finished in THIS cycle
            if (rs.cyclesRemaining == 0) {
                // Finish execution in THIS cycle
                int result = executeInstruction(rs);
                ROBEntry& entry = reorderBuffer[rs.robTag];
                entry.value = to16Bit(result);
                // CRITICAL: Mark as execution finished but DON'T set ready=true yet
                entry.endExecCycle = currentCycle;
                instructionTimings[entry.instructionAddress].endExecCycle = currentCycle;
                
                // Move to WRITING state - will set ready=true in NEXT cycle
                rs.state = WRITING;
                rs.cyclesRemaining = 1; // Ensure it's 0
                
                if (debugMode) {
                    cout << "Cycle " << currentCycle << ": Finish exec RS " << rs.id 
                         << ", result=" << entry.value << " (will be ready next cycle)\n";
                }
            }
        }
    }
}
int executeInstruction(ReservationStation& rs) {
        int result = 0;
        
        switch (rs.type) {
            case LOAD: {
                int address = to16Bit(rs.vj + rs.offset);
                // result = memory[address];
                if (address < 0 || address >= MEMORY_SIZE) {
    cerr << "Memory access out of bounds at address " << address << "\n";
    result = 0;
} else {
    result = memory[address];
}
                if (debugMode) {
                    cout << "    LOAD: " << rs.vj << " + " << rs.offset << " = addr " 
                         << address << " -> " << result << "\n";
                }
                break;
            }
            case STORE: {
                int address = to16Bit(rs.vj + rs.offset);
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
                    if (debugMode) cout << "    BEQ misprediction!\n";
                }
                break;
            }
            case CALL: {
                // Store PC+1 in R1
                result = reorderBuffer[rs.robTag].instructionAddress + 1;
                // Calculate target address
                int target = result + rs.offset;
                reorderBuffer[rs.robTag].value = to16Bit(target);  // Store target in value field
                if (debugMode) {
                    cout << "    CALL: return addr=" << result << ", target=" << target << "\n";
                }
                break;
            }
            case RET: {
                       if (rs.vj == 0) {
                // Check if R1 has a valid value
                if (registerFile[1].valid) {
                    result = registerFile[1].value;
                    cout << "RET FIX: Using R1 value directly: " << result << endl;
                } else {
                    result = rs.vj;
                }
            } else {
                result = rs.vj;
            }
            
            if (debugMode) {
                cout << "    RET: jumping to " << result << "\n";
            }
            break;
        }
            case ADD:
                result = rs.vj + rs.vk;
                if (debugMode) cout << "    ADD: " << rs.vj << " + " << rs.vk << " = " << result << "\n";
                break;
            case SUB:
                result = rs.vj - rs.vk;
                if (debugMode) cout << "    SUB: " << rs.vj << " - " << rs.vk << " = " << result << "\n";
                break;
            case NAND:
                result = ~(rs.vj & rs.vk);
                if (debugMode) cout << "    NAND: " << rs.vj << " NAND " << rs.vk << " = " << result << "\n";
                break;
            case MUL: {
                int32_t product = (int32_t)rs.vj * (int32_t)rs.vk;
                result = product & 0xFFFF;
                if (debugMode) cout << "    MUL: " << rs.vj << " * " << rs.vk << " = " << result << "\n";
                break;
            }
            default:
                break;
        }
        
        return result;
    }
    
    int getCyclesForType(InstructionType type) {
        switch (type) {
            case LOAD: return 6;      // 2 (address) + 4 (memory)
            case STORE: return 2;     // 2 (address) + 4 (memory)
            case BEQ: return 1;       // Compare + compute target
            case CALL: return 1;      // Compute target + store return
            case RET: return 1;       // Branch to R1
            case ADD: return 2;
            case SUB: return 2;
            case NAND: return 1;
            case MUL: return 12;
            default: return 1;
        }
    }
void writeResultStage() {
    for (auto& rs : reservationStations) {
        if (rs.state == WRITING) {
            // Only process if delay is already finished (was 0 before this cycle)
            if (rs.cyclesRemaining == 0) {
                ROBEntry& robEntry = reorderBuffer[rs.robTag];
                robEntry.ready = true;
                
                // Check for speculative flush
                if (robEntry.isSpeculative && robEntry.speculativeBranchTag != -1) {
                    ROBEntry& branchEntry = reorderBuffer[robEntry.speculativeBranchTag];
                    if (branchEntry.type == BEQ && branchEntry.ready && branchEntry.value == 1) {
                        if (debugMode) {
                            cout << "Cycle " << currentCycle << ": SKIP write RS " << rs.id 
                                 << " ROB " << rs.robTag << " - speculative on mispredicted branch\n";
                        }
                        rs.state = IDLE;
                        rs.robTag = -1;
                        rs.startedExecution = false;
                        continue;
                    }
                }
                
                // Broadcast and write
                broadcastResult(rs.robTag, robEntry.value);
                robEntry.writeCycle = currentCycle;
                instructionTimings[robEntry.instructionAddress].writeCycle = currentCycle;
                
                if (debugMode) {
                    cout << "Cycle " << currentCycle << ": Write result RS " << rs.id 
                         << " ROB " << rs.robTag << " = " << robEntry.value 
                         << " (exec finished cycle " << robEntry.endExecCycle << ")\n";
                }
                
                // Clear RS
                rs.state = IDLE;
                rs.robTag = -1;
                rs.startedExecution = false;
            } else {
                // Still waiting - decrement counter
                rs.cyclesRemaining--;
            }
        }
    }
}

void broadcastResult(int robTag, int value) {
    // Update reservation stations
    for (auto& rs : reservationStations) {
        if (rs.qj == robTag) {
            rs.vj = value;
            rs.qj = -1;
            if (debugMode) {
                cout << "    Broadcast to RS" << rs.id << " (type: " 
                     << instructionTypeToString(rs.type) << "): vj = " << value << "\n";
            }
        }
        if (rs.qk == robTag) {
            rs.vk = value;
            rs.qk = -1;
            if (debugMode) {
                cout << "    Broadcast to RS" << rs.id << " (type: " 
                     << instructionTypeToString(rs.type) << "): vk = " << value << "\n";
            }
        }
    }
    
    // ALSO: If any instruction is waiting for a ROB that has been freed
    // (happens when ROB entry gets reused), check current register values
    for (auto& rs : reservationStations) {
        if (rs.qj != -1) {
            // Check if this ROB tag still exists
            bool exists = false;
            for (const auto& entry : reorderBuffer) {
                if (entry.tag == rs.qj && entry.type != INVALID) {
                    exists = true;
                    break;
                }
            }
            
            if (!exists) {
                // ROB was freed, check register directly
                // RET reads from R1, others might read from different regs
                if (rs.type == RET) {
                    if (registerFile[1].valid) {
                        rs.vj = registerFile[1].value;
                        rs.qj = -1;
                        cout << "FIXED: RET now has vj=" << rs.vj << endl;
                    }
                }
            }
        }
    }
}
void simulate() {
    nextInstructionToIssue = 0;
    branchMispredicted = false;
    branchTarget = -1;
    int maxInstructions = instructionQueue.size();
       // ADD THIS DEBUG
    cout << "\n=== DEBUG: Instruction Queue Contents ===\n";
    for (int i = 0; i < maxInstructions; i++) {
        const Instruction& instr = instructionQueue[i];
        cout << "Index " << i << ": address=" << instr.address 
             << ", type=" << instructionTypeToString(instr.type)
             << ", label='" << instr.label << "'\n";
        cout << "  rA=" << instr.rA << ", rB=" << instr.rB 
             << ", rC=" << instr.rC << ", offset=" << instr.offset << "\n";
        cout << "  Valid: " << (instr.type != INVALID ? "YES" : "NO") << "\n";
        cout << "---\n";
    }
    cout << "=========================================\n\n";
    cout << "\n=== Starting Simulation ===\n";
    cout << "Total instructions: " << maxInstructions << "\n";
    cout << "Using ALWAYS-NOT-TAKEN branch predictor\n";
    cout << "ROB size: " << ROB_SIZE << "\n";
    cout << "Starting address: " << programStartAddress << "\n\n";

    bool allWorkDone = false;
    int stallCounter = 0;
    int lastCommitCount = 0;
    
    while (!allWorkDone) {
        currentCycle++;
        totalCycles = currentCycle;

        if (debugMode) {
            cout << "\n--- Cycle " << currentCycle << " ---\n";
        }
        
        // Process stages in reverse order (commit -> execute -> write -> issue)
        commitStage();           // Commit first
        executeStage();          // Then execute
        writeResultStage();      // Then write results
        
        // Issue next instruction
        if (nextInstructionToIssue < maxInstructions && !branchMispredicted) {
            if (issueInstruction(instructionQueue[nextInstructionToIssue])) {
                nextInstructionToIssue++;
            }
        }

        // Handle branch misprediction redirect (AFTER all pipeline stages)
        if (branchMispredicted) {
            if (debugMode) {
                cout << "  Redirecting fetch to address " << branchTarget << "\n";
            }
            
            // Flush the instruction queue and restart from target
            nextInstructionToIssue = 0;
            
            // Find the instruction at the target address
            bool found = false;
            for (size_t i = 0; i < instructionQueue.size(); i++) {
                if (instructionQueue[i].address == branchTarget) {
                    nextInstructionToIssue = i;
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                // If target not found in current instructions, we're done
                nextInstructionToIssue = maxInstructions;
            }
            
            branchMispredicted = false;
        }

        if (debugMode) {
            printReservationStations();
            printRegisterStatus(); 
            printROB();
            cout << string(50, '-') << "\n";  // Separator
        }

        // ====== FIXED TERMINATION LOGIC ======
        // Check if all instructions are done
        if (instructionsCompleted >= maxInstructions ) {
            allWorkDone = true;
            if (debugMode) {
                cout << "TERMINATING: All " << maxInstructions << " instructions completed.\n";
            }
            break;  // Exit the loop immediately
        }
        
        // Check if we should stop (all instructions issued AND completed)
        bool robEmpty = true;
        for (const auto& entry : reorderBuffer) {
            if (entry.type != INVALID) {
                robEmpty = false;
                break;
            }
        }
        
        bool rsIdle = true;
        for (const auto& rs : reservationStations) {
            if (rs.state != IDLE) {
                rsIdle = false;
                break;
            }
        }
        
        // If all instructions issued AND ROB empty AND RS idle, we're done
        if (nextInstructionToIssue >= maxInstructions && robEmpty && rsIdle) {
            allWorkDone = true;
            if (debugMode) {
                cout << "TERMINATING: All instructions issued, ROB empty, RS idle.\n";
            }
            break;
        }

        // Stall detection - only if we haven't finished all instructions
        if (instructionsCompleted == lastCommitCount) {
            stallCounter++;
            if (stallCounter > 50) {
                cout << "\nWarning: Detected stall for 50 cycles. Possible deadlock.\n";
                cout << "ROB Head: " << robHead << ", Tail: " << robTail << "\n";
                cout << "Instructions completed: " << instructionsCompleted << "/" << maxInstructions << "\n";
                cout << "Next instruction to issue: " << nextInstructionToIssue << "/" << maxInstructions << "\n";
                
                // Print ROB state for debugging
                bool anyValid = false;
                for (int i = 0; i < ROB_SIZE; i++) {
                    const ROBEntry& e = reorderBuffer[i];
                    if (e.type != INVALID) {
                        anyValid = true;
                        cout << "  ROB[" << i << "]: " << instructionTypeToString(e.type)
                             << ", ready=" << (e.ready ? "Y" : "N")
                             << ", dest=R" << e.destination << "\n";
                    }
                }
                if (!anyValid) {
                    cout << "  ROB is completely empty!\n";
                }
                
                // Force termination if stall detected
                allWorkDone = true;
                break;
            }
        } else {
            stallCounter = 0;
            lastCommitCount = instructionsCompleted;
        }

        // Safety limit
        if (currentCycle > 100) {
            cout << "\nWarning: Stopping after 100 cycles\n";
            cout << "Instructions completed: " << instructionsCompleted << "/" << maxInstructions << "\n";
            allWorkDone = true;
            break;
        }
    }
    // Issue next instruction
if (nextInstructionToIssue < maxInstructions && !branchMispredicted) {
    if (debugMode && currentCycle > 20) {
        cout << "  Trying to issue instruction " << nextInstructionToIssue 
             << " at address " << instructionQueue[nextInstructionToIssue].address
             << ": " << instructionQueue[nextInstructionToIssue].label << "\n";
    }
    
    if (issueInstruction(instructionQueue[nextInstructionToIssue])) {
        nextInstructionToIssue++;
    } else {
        if (debugMode && currentCycle > 20) {
            cout << "  Failed to issue instruction " << nextInstructionToIssue << "\n";
        }
    }
}

if (debugMode && currentCycle > 20) {
    cout << "  nextInstructionToIssue = " << nextInstructionToIssue 
         << "/" << maxInstructions << "\n";
    cout << "  branchMispredicted = " << (branchMispredicted ? "true" : "false") << "\n";
}
    cout << "\nSimulation completed in " << totalCycles << " cycles.\n";
}
void commitStage() {
    // If we're in the middle of committing a store, handle that first
    if (commitCyclesRemaining > 0) {
        commitCyclesRemaining--;
        if (commitCyclesRemaining == 0 && committingROBIndex != -1) {
            // Finish the store commit
            ROBEntry& entry = reorderBuffer[committingROBIndex];
            
            // Write to memory
            if (entry.destination >= 0 && entry.destination < (int)memory.size()) {
                memory[entry.destination] = entry.value;
                if (debugMode) {
                    cout << "Cycle " << currentCycle << ": Store to memory[" 
                         << entry.destination << "] = " << entry.value << "\n";
                }
            }
            
            // Mark committed
            entry.commitCycle = currentCycle;
            if (entry.instructionAddress >= 0 && entry.instructionAddress < (int)instructionTimings.size()) {
                instructionTimings[entry.instructionAddress].commitCycle = currentCycle;
                instructionTimings[entry.instructionAddress].committed = true;
            }
            instructionsCompleted++;
            
            // Free ROB entry
            entry.state = ISSUED;
            entry.type = INVALID;
            entry.ready = false;
            entry.destination = -1;
            entry.isSpeculative = false;
            entry.speculativeBranchTag = -1;
            
            // Advance head pointer to next entry
            robHead = (committingROBIndex + 1) % ROB_SIZE;
            
            committingROBIndex = -1;
        }
        return;
    }
    
    // Try to commit the instruction at the HEAD of ROB
    ROBEntry& entry = reorderBuffer[robHead];
    
    // FIXED: Only check if entry is valid and ready
    if (entry.type != INVALID && entry.ready) {
        // Check for branch misprediction FIRST
        if (entry.type == BEQ && entry.value == 1) {
            if (debugMode) {
                cout << "Cycle " << currentCycle << ": Commit ROB " << robHead 
                     << " BEQ MISPREDICTION - FLUSHING\n";
            }
            branchMispredictions++;
            
            // Flush all speculative instructions
            flushSpeculativeInstructions(robHead);
            
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
                cout << "Branch taken: PC=" << branchPC << "+1+offset=" 
                     << branchOffset << " -> target address " << targetAddress << "\n";
            }
            
            // Set redirection flags for next cycle
            branchMispredicted = true;
            branchTarget = targetAddress;
            
            // Commit the branch itself
            entry.commitCycle = currentCycle;
            if (entry.instructionAddress >= 0 && entry.instructionAddress < (int)instructionTimings.size()) {
                instructionTimings[entry.instructionAddress].commitCycle = currentCycle;
                instructionTimings[entry.instructionAddress].committed = true;
            }
            instructionsCompleted++;
            
            // Clear unresolved branch tracking
            if (lastUnresolvedBranchROB == robHead) {
                lastUnresolvedBranchROB = -1;
            }
            
            // Free ROB entry and advance head
            entry.state = ISSUED;
            entry.type = INVALID;
            entry.ready = false;
            entry.isSpeculative = false;
            entry.speculativeBranchTag = -1;
            
            robHead = (robHead + 1) % ROB_SIZE;
            
            // Check if we need to wrap around to find next valid entry
            if (reorderBuffer[robHead].type == INVALID) {
                int nextValid = findNextValidROBEntry(robHead);
                if (nextValid != -1) {
                    robHead = nextValid;
                } else {
                    // ROB is empty
                    robHead = 0;
                    robTail = 0;
                }
            }
            
            return;
        }
        
        // Check for CALL/RET redirection
        if (entry.type == CALL || entry.type == RET) {
            if (debugMode) {
                cout << "Cycle " << currentCycle << ": Commit ROB " << robHead << " " 
                     << instructionTypeToString(entry.type) << " - REDIRECTING\n";
            }
            
            // Calculate target address
            int targetAddress;
            if (entry.type == CALL) {
                targetAddress = entry.instructionAddress + 1 + entry.callOffset;
            } else {
                targetAddress = entry.value;
            }
            
            // Set redirection
            branchMispredicted = true;
            branchTarget = targetAddress;
            
            if (debugMode) {
                cout << "  Redirecting to address " << targetAddress << "\n";
            }
        }
        
        // Skip committing if speculative on mispredicted branch
        if (entry.isSpeculative && entry.speculativeBranchTag != -1) {
            ROBEntry& branchEntry = reorderBuffer[entry.speculativeBranchTag];
            if (branchEntry.type == BEQ && branchEntry.ready && branchEntry.value == 1) {
                if (debugMode) {
                    cout << "Cycle " << currentCycle << ": SKIP commit ROB " 
                         << robHead << " - speculative on mispredicted branch\n";
                }
                
                // Still need to free the entry
                entry.state = ISSUED;
                entry.type = INVALID;
                entry.ready = false;
                entry.isSpeculative = false;
                entry.speculativeBranchTag = -1;
                
                robHead = (robHead + 1) % ROB_SIZE;
                
                // Find next valid entry if current is invalid
                if (robHead != robTail && reorderBuffer[robHead].type == INVALID) {
                    int nextValid = findNextValidROBEntry(robHead);
                    if (nextValid != -1) {
                        robHead = nextValid;
                    }
                }
                return;
            }
        }
        
        // Handle STORE - takes 4 cycles to commit
        if (entry.type == STORE) {
            if (debugMode) {
                cout << "Cycle " << currentCycle << ": Start commit ROB " << robHead << " " 
                     << instructionTypeToString(entry.type) << " (4 cycles for memory write)\n";
            }
            commitCyclesRemaining = 4;
            committingROBIndex = robHead;
            return;
        }
        
        // Normal commit (instant for non-STORE instructions)
        if (debugMode) {
            cout << "Cycle " << currentCycle << ": Commit ROB " << robHead << " " 
                 << instructionTypeToString(entry.type) << "\n";
        }
        
        // Update architectural state at commit time
        if (entry.destination >= 0 && entry.destination < NUM_REGISTERS) {
            if (registerFile[entry.destination].robTag == robHead) {
                registerFile[entry.destination].value = entry.value;
                registerFile[entry.destination].valid = true;
                registerFile[entry.destination].robTag = -1;
                if (debugMode) {
                    cout << "  Update R" << entry.destination << " = " << entry.value << "\n";
                }
            }
        }
        
        // Mark committed and update timing
        entry.commitCycle = currentCycle;
        if (entry.instructionAddress >= 0 && entry.instructionAddress < (int)instructionTimings.size()) {
            instructionTimings[entry.instructionAddress].commitCycle = currentCycle;
            instructionTimings[entry.instructionAddress].committed = true;
        }
        instructionsCompleted++;
        
        // Clear unresolved branch tracking if this was a branch
        if (entry.type == BEQ && lastUnresolvedBranchROB == robHead) {
            lastUnresolvedBranchROB = -1;
        }
        
        // Free ROB entry
        entry.state = ISSUED;
        entry.type = INVALID;
        entry.ready = false;
        entry.destination = -1;
        entry.isSpeculative = false;
        entry.speculativeBranchTag = -1;
        
        // Advance head pointer
        robHead = (robHead + 1) % ROB_SIZE;
        
        // Skip any invalid entries at the new head
        while (robHead != robTail && reorderBuffer[robHead].type == INVALID) {
            robHead = (robHead + 1) % ROB_SIZE;
        }
        
        // If ROB is empty, reset pointers
        if (robHead == robTail && reorderBuffer[robHead].type == INVALID) {
            robHead = 0;
            robTail = 0;
        }
    }
    // Handle case where head entry is invalid but there might be valid entries
    else if (entry.type == INVALID && robHead != robTail) {
        // Try to find next valid entry
        int nextValid = findNextValidROBEntry(robHead);
        if (nextValid != -1) {
            robHead = nextValid;
        } else {
            // All entries are invalid, reset
            robHead = 0;
            robTail = 0;
        }
    }
}

// Helper function to find next valid ROB entry
int findNextValidROBEntry(int start) {
    int current = (start + 1) % ROB_SIZE;
    while (current != start) {
        if (reorderBuffer[current].type != INVALID) {
            return current;
        }
        current = (current + 1) % ROB_SIZE;
    }
    return -1;  // No valid entries found
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
    // Add these after line with instructionTypeToString function
string stateToString(ReservationStationState state) {
    switch(state) { 
        case IDLE: return "IDLE"; 
        case BUSY: return "BUSY"; 
        case EXECUTING: return "EXEC"; 
        case WRITING: return "WRITE"; 
        default: return "?"; 
    }
}

string robStateToString(ROBState state) {
    switch(state) { 
        case ISSUED: return "ISSUED"; 
        case EXECUTING_ROB: return "EXEC";  // <- Fixed: EXECUTING_ROB not EXECUTINGROB
        case WRITTEN: return "WRITTEN"; 
        case COMMITTED: return "COMMIT"; 
        default: return "?"; 
    }
}


    void loadDefaultProgram() {
        cout << "Loading default test program...\n";
        
        instructionQueue.clear();
        instructionQueue.push_back(parseInstruction("LOAD R1, 0(R0)", 0));
        instructionQueue.push_back(parseInstruction("ADD R2, R1, R1", 1));
        instructionQueue.push_back(parseInstruction("STORE R2, 4(R0)", 2));
        instructionQueue.push_back(parseInstruction("BEQ R1, R0, 2", 3));
        instructionQueue.push_back(parseInstruction("ADD R3, R2, R1", 4));
        instructionQueue.push_back(parseInstruction("STORE R3, 8(R0)", 5));
        
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
        
        // Get starting address
        cout << "Enter program starting address (0-" << MEMORY_SIZE-1 << "): ";
        cin >> programStartAddress;
        programStartAddress = programStartAddress & 0xFFFF;
        
        instructionQueue.clear();
        
        // Load instructions
        ifstream progFile("program.txt");
        if (!progFile.is_open()) {
            cout << "Error: Cannot open program.txt\n";
            cout << "Please create a program.txt file with instructions.\n";
            cout << "Using default program instead.\n";
            loadDefaultProgram();
            return;
        }
        
        string line;
        int address = programStartAddress;
        while (getline(progFile, line)) {
            // Skip empty lines
            bool empty = true;
            for (char c : line) {
                if (!isspace(static_cast<unsigned char>(c))) { 
                    empty = false; 
                    break; 
                }
            }
            if (!empty && line[0] != '#') {  // Skip comments
                instructionQueue.push_back(parseInstruction(line, address++));
            }
        }
        progFile.close();
        
        if (instructionQueue.empty()) {
            cout << "No instructions found in program.txt\n";
            cout << "Using default program instead.\n";
            loadDefaultProgram();
            return;
        }
        
        instructionTimings.clear();
        instructionTimings.resize(MEMORY_SIZE);
      for (int i = 0; i < instructionQueue.size(); i++) {
    instructionTimings[i].address = instructionQueue[i].address;
}
        // for (int i = 0; i < (int)instructionQueue.size(); i++) {
        //     instructionTimings[instructionQueue[i].address].address = instructionQueue[i].address;
        // }
        
        // Clear and load memory
        for (int i = 0; i < MEMORY_SIZE; i++) {
            memory[i] = 0;
        }
        
        ifstream memFile("memory.txt");
        if (memFile.is_open()) {
            int addr, value;
            while (memFile >> addr >> value) {
                if (addr >= 0 && addr < MEMORY_SIZE) {
                    memory[addr] = value & 0xFFFF;
                }
            }
            memFile.close();
        } else {
            cout << "Warning: memory.txt not found (starting with all zeros)\n";
        }
        
        // Reset registers
        for (int i = 0; i < NUM_REGISTERS; i++) {
            registerFile[i].value = 0;
            registerFile[i].valid = true;
            registerFile[i].robTag = -1;
        }
        registerFile[0].value = 0;
        
        cout << "Program loaded from files:\n";
        cout << "Starting address: " << programStartAddress << "\n";
        for (int i = 0; i < (int)instructionQueue.size(); i++) {
            cout << instructionQueue[i].address << ": " << instructionQueue[i].label << "\n";
        }
        cout << "\n=== PARSING DEBUG ===\n";
while (getline(progFile, line)) {
    // Skip empty lines
    bool empty = true;
    for (char c : line) {
        if (!isspace(c)) { 
            empty = false; 
            break; 
        }
    }
    
    if (!empty && line[0] != '#') {
        cout << "Line: '" << line << "'\n";
        Instruction instr = parseInstruction(line, address);
        cout << "  Type: " << instructionTypeToString(instr.type) 
             << ", Label: '" << instr.label << "'\n";
        cout << "  rA: " << instr.rA << ", rB: " << instr.rB 
             << ", offset: " << instr.offset << "\n";
        
        if (instr.type != INVALID) {
            instructionQueue.push_back(instr);
            address++;
            cout << "  -> ADDED to queue\n";
        } else {
            cout << "  -> SKIPPED (invalid)\n";
        }
        cout << "---\n";
    }
}
cout << "====================\n\n";
    }
    
    void printResults() {
        cout << "\n=== SIMULATION RESULTS ===\n\n";
        
        cout << "Instruction Timeline:\n";
        cout << "============================================================================\n";
        cout << left << setw(5) << "Addr" << setw(20) << "Instruction" 
             << setw(8) << "Issue" << setw(8) << "Exec" 
             << setw(8) << "Finish" << setw(8) << "Write" 
             << setw(8) << "Commit" << "\n";
        cout << "============================================================================\n";
          cout << left << setw(5) << "Addr" << setw(20) << "Instruction" 
         << setw(8) << "Issue" << setw(8) << "Exec" 
         << setw(8) << "Finish" << setw(8) << "Write" 
         << setw(8) << "Commit" << "\n";
    cout << "============================================================================\n";
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
        cout << "(Using ALWAYS-NOT-TAKEN branch predictor)\n";
        
        cout << "\nRegister File State:\n";
        for (int i = 0; i < NUM_REGISTERS; i++) {
            cout << "R" << i << ": " << registerFile[i].value;
            if (!registerFile[i].valid) {
                cout << " (Waiting for ROB " << registerFile[i].robTag << ")";
            }
            cout << "\n";
        }
        
        cout << "\nMemory Contents (first 10 locations):\n";
        for (int i = 0; i < 10; i++) {
            cout << "M[" << i << "]: " << memory[i] << "\n";
        }
    }
    // Add these 3 functions right after printResults() function ends
void printReservationStations() {
    cout << "\n=== RESERVATION STATIONS (Cycle " << currentCycle << ") ===\n";
    cout << left << setw(4) << "ID" << setw(8) << "Type" << setw(8) << "State" 
         << setw(4) << "A" << setw(4) << "B" << setw(6) << "Vj" << setw(6) << "Vk" 
         << setw(6) << "Qj" << setw(6) << "Qk" << setw(4) << "Rem" << endl;
    for (const auto& rs : reservationStations) {
        if (rs.state != IDLE) {  // Only show busy stations
            cout << left << setw(4) << rs.id << setw(8) << instructionTypeToString(rs.type)
                 << setw(8) << stateToString(rs.state) << setw(4) << rs.rA 
                 << setw(4) << rs.rB << setw(6) << rs.vj << setw(6) << rs.vk 
                 << setw(6) << (rs.qj != -1 ? to_string(rs.qj) : "-") 
                 << setw(6) << (rs.qk != -1 ? to_string(rs.qk) : "-")
                 << setw(4) << rs.cyclesRemaining << endl;
        }
    }
}

void printRegisterStatus() {
    cout << "\n=== REGISTER STATUS ===\n";
    cout << left << setw(2) << "R" << setw(6) << "Value" << setw(6) << "Valid" << setw(4) << "ROB" << endl;
    for (int i = 0; i < registerFile.size(); ++i) {  // <- Fixed: use registerFile.size()
        cout << left << setw(2) << i << setw(6) << registerFile[i].value 
             << setw(6) << (registerFile[i].valid ? "Y" : "N") 
             << setw(4) << registerFile[i].robTag << endl;
    }
}


void printROB() {
    cout << "\n=== REORDER BUFFER ===\n";
    cout << left << setw(3) << "Tag" << setw(10) << "State" << setw(6) << "Type" 
         << setw(5) << "Dest" << setw(5) << "Value" << setw(4) << "Rdy" << endl;
    for (const auto& entry : reorderBuffer) {
        if (entry.type != INVALID) {
            cout << left << setw(3) << entry.tag << setw(10) << robStateToString(entry.state)
                 << setw(6) << instructionTypeToString(entry.type) << setw(5) << entry.destination 
                 << setw(5) << entry.value << setw(4) << (entry.ready ? "Y" : "N") << endl;
        }
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
    cout << "==============================================\n";
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

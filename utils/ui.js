// user interaction
let showMethodButtons = [];
let rigMethodButtons = [];
let showMethods = [
  { name: "Rig", value: SoftBodyShowMethod.RIG },
  { name: "Surface", value: SoftBodyShowMethod.SURFACE }
];
let rigMethods = [
  { name: "Perimeter", value: SoftBodyRigMethod.PERIMETER },
  { name: "Grid", value: SoftBodyRigMethod.GRID },
  // { name: "Half grid", value: SoftBodyRigMethod.HALFGRID },
];


let selectedRigMethod = SoftBodyRigMethod.GRID;
let selectedShowMethod = SoftBodyShowMethod.RIG;


function createSoftBodyUI() {
    const uiParent = select('#ui-container');
  
    // Show method buttons
    let showDiv = createDiv('Show Method:');
    showDiv.addClass('centered-div');
    showDiv.parent(uiParent);
    showMethodButtons = [];
    showMethods.forEach((method, idx) => {
      let btn = createButton(method.name);
      btn.parent(showDiv);
      btn.mousePressed(() => {
        selectedShowMethod = method.value;
        showMethodButtons.forEach(b => b.removeClass('selected'));
        btn.addClass('selected');
      });
      showMethodButtons.push(btn);
      if (idx === 0) btn.addClass('selected');
    });
  
    // Rig method buttons
    let rigDiv = createDiv('Rig Method:');
    rigDiv.addClass('centered-div');
    rigDiv.parent(uiParent);
    rigMethodButtons = [];
    rigMethods.forEach((method, idx) => {
      let btn = createButton(method.name);
      btn.parent(rigDiv);
      btn.mousePressed(() => {
        selectedRigMethod = method.value;
        rigMethodButtons.forEach(b => b.removeClass('selected'));
        btn.addClass('selected');
        restartSimulation();
      });
      rigMethodButtons.push(btn);
      if (idx === 0) btn.addClass('selected');
    });
  
    // Restart button
    let restartBtn = createButton('Restart Simulation');
    restartBtn.addClass('centered-div');
    restartBtn.parent(uiParent);
    restartBtn.mousePressed(() => {
      restartSimulation();
    });
  }
  
  function restartSimulation() {
    softBody = new Soft_body(350, 10, SOFTBODY_WIDTH, SOFTBODY_HEIGHT, SPRING_CONSTRAINT_STIFFNESS, SOLVER_ITERATIONS, selectedRigMethod, selectedShowMethod);
    hashGrid = new Grid(2 * SOFTBODY_SIZE, softBody.masses.length);
    // If you need to reset other state (like hashGrid), do it here as well
  }
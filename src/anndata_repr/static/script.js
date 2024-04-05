let root = document.getElementById("__ID__");

let sections = {};
for (let section of root.querySelectorAll(".ad-section-item")) {
  let input = section.querySelector("input[type=checkbox]");
  let name = input?.dataset?.anndata;
  if (!name) continue;
  sections[name] = input;
}

let svg = root.querySelector(".ad-svg");
for (let grp of svg.querySelectorAll("g")) {
  let id = grp?.firstChild?.id;
  let input = sections[id === "X" ? "layers" : id];
  if (!input) continue;
  let label = input.nextElementSibling;
  let fill = grp.querySelector("rect").getAttribute("fill");

  grp.addEventListener("mouseenter", () => {
    label.style.color = fill;
  });

  grp.addEventListener("mouseleave", () => {
    label.style.color = "";
  });

  grp.addEventListener("click", () => {
    input.checked = !input.checked;
  });
}

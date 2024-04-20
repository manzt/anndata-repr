class Section {
  #input;
  #label;
  constructor(root, name, rect) {
    this.root = root;
    this.name = name;
    this.rect = rect;
    this.#input = root.querySelector("input[type=checkbox]");
    this.#label = root.querySelector("label");
    this.#setupMouseEvents();
  }
  expand() {
    this.#input.checked = true;
  }
  collapse() {
    this.#input.checked = false;
  }
  #setupMouseEvents() {
    if (!this.rect) return;
    // Add click event to the label
    const g = this.rect.parentElement;
    const fill = this.rect.getAttribute("fill");
    g.addEventListener("mouseenter", () => {
      if (fill === "transparent") return;
      this.#label.style.color = fill;
    });
    g.addEventListener("mouseleave", () => {
      this.#label.style.color = "";
    });
    g.addEventListener("click", () => {
      if (!this.#input.disabled) {
        this.#input.checked = !this.#input.checked;
      }
    });
  }
  varItems() {
    return [...this.root.querySelectorAll(".ad-var-item")].map((el) => ({
      hide: () => el.style.display = "none",
      show: () => el.style.display = "",
      textContent: el.querySelector(".ad-var-name").textContent,
    }));
  }
}

function gatherSections(root) {
  const svg = root.querySelector(".ad-svg");
  /** @type {Array<Section>} */
  const sections = [];
  for (const section of root.querySelectorAll(".ad-section-item")) {
    const input = section.querySelector("input[type=checkbox]");
    const name = input?.dataset?.anndata;
    if (!name) continue;
    const rect = svg.querySelector(`#${name === "layers" ? "X" : name}`);
    sections.push(new Section(section, name, rect));
  }
  return sections;
}

// Simple matching function, case-insensitive.
// Returns true if `search` is a prefix of `content` or any part of `content`
// separated by "_"
function matches(content, search) {
  content = content.toLowerCase();
  return content.startsWith(search) ||
    content.split("_").some((part) => part.startsWith(search));
}

function addSearch(sections) {
  const form = root.querySelector(".ad-search");
  const input = form.querySelector("input");

  // Focus entire text input on focus event
  input.addEventListener("focus", () => {
    input.select();
  });

  input.addEventListener("input", () => {
    const search = input.value.trim().toLowerCase();

    // Show/hide variables
    if (!search) {
      for (const section of sections) {
        section.collapse();
        for (const item of section.varItems()) {
          item.show();
        }
      }
      return;
    }

    for (const section of sections) {
      section.collapse();
      for (const item of section.varItems()) {
        if (matches(item.textContent, search)) {
          item.show();
          section.expand();
        } else {
          item.hide();
        }
      }
    }
  });

  form.addEventListener("submit", (e) => {
    e.preventDefault();
    input.blur();
  });
}

const root = document.getElementById("__ID__");
const sections = gatherSections(root);
addSearch(sections);

function main() {
    const root = document.getElementById("__ID__")
    console.log(root)
    const svg = root.querySelector(".adata-overview")
    console.log(svg)
    if (!svg) return

    const blocks = svg.querySelectorAll(".ad-block")
    const elements = root.querySelectorAll(".ad-section-summary-in")

    for (let b of blocks) {
        let anndataName = b.getAttribute("id")
        b.addEventListener("mouseover", () => {
            change_opacity(blocks, b, 0.5)
            addNameBackground(b)
        })
        b.addEventListener("mouseout", () => {
            change_opacity(blocks, 0, 1)
            removeNameBackground(b)
        })
        b.addEventListener("click", (event) => {
            for (let e of elements) {
                if (e.dataset?.anndata !== anndataName) {
                    e.checked = false
                } else {
                    e.checked = true
                }
            }
        })
    }

    function removeNameBackground(block) {
        block.querySelectorAll(".text-background").forEach((rect) => {
            // Remove the rect when not hovered
            block.removeChild(rect)
        })
    }

    function addNameBackground(blockElement) {
        blockElement.querySelectorAll("text").forEach((text) => {
            let bbox = text.getBBox()
            let padding = 2 // Adjust padding around the text
            let rect = document.createElementNS("http://www.w3.org/2000/svg", "rect")
            rect.setAttribute("x", bbox.x - padding)
            rect.setAttribute("y", bbox.y - padding)
            rect.setAttribute("width", bbox.width + 2 * padding)
            rect.setAttribute("height", bbox.height + 2 * padding)
            rect.setAttribute("class", "text-background")
            rect.setAttribute("pointer-events", "none")

            // Insert the rect before the text node in the SVG.
            blockElement.insertBefore(rect, text)
        })
    }

    function change_opacity(blocks, except_block, opacity) {
        for (let b of blocks) {
            if (b !== except_block) {
                b.style.opacity = opacity
                b.style.transitionDuration = "250ms"
            }
        }
    }
}
main()

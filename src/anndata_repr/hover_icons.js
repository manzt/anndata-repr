const svg_id = "__REPLACE_ME__"

const svg = document.getElementById(svg_id);

const blocks = svg.querySelectorAll(".block")


for (let b of blocks) {
    b.addEventListener("mouseover", () => {
        change_opacity(blocks, b, 0.5)
    })
    b.addEventListener("mouseout", () => {
        change_opacity(blocks, 0, 1)
    })
    b.addEventListener("click", () => {

        const elements = document.getElementsByClassName('ad-section-summary-in')
        for (let e of elements) {
            if(e.dataset.blockname !== b.getAttribute("id")){
                e.checked = false
            } else {
                e.checked = true
            }
        }

        
    })
}

function change_opacity(blocks, except_block, opacity) {
    for (let b of blocks) {
        if (b !== except_block) {
            b.style.opacity = opacity;
            b.style.transitionDuration = "250ms";
        }
    }
}

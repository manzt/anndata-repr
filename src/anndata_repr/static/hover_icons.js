function main() {
    const root = document.getElementById("__ID__");
    const svg = root.querySelector(".adata-overview");
    if (!svg) return;

    const blocks = svg.querySelectorAll(".ad-block")
    const elements = root.querySelectorAll(".ad-section-summary-in")

    for (let b of blocks) {
        let anndataName = b.getAttribute("id");
        console.log(anndataName);

        b.addEventListener("mouseover", () => {
            change_opacity(blocks, b, 0.5)
        })
        b.addEventListener("mouseout", () => {
            change_opacity(blocks, 0, 1)
        })
        b.addEventListener("click", (event) => {
            for (let e of elements) {
                if (e.dataset?.anndata !== b.getAttribute("id")){
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
}
main();

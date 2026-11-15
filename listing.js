// Decorates blog/tidbit listings: entry numbers, size pips and reading time.
document.addEventListener("DOMContentLoaded", function () {
    const sizes = { nibble: 1, tidbit: 2, chunk: 3 };
    const posts = Array.from(document.querySelectorAll(".quarto-post"));

    // number entries chronologically (oldest = 01), independent of display order
    const byDate = posts.slice().sort((a, b) =>
        Number(a.dataset.listingDateSort) - Number(b.dataset.listingDateSort));

    byDate.forEach((post, i) => {
        const meta = post.querySelector(".metadata a") || post.querySelector(".metadata");
        if (!meta) return;

        const number = document.createElement("div");
        number.className = "listing-number";
        number.innerHTML = `<span class="listing-number-hash">#</span>${String(i + 1).padStart(2, "0")}`;
        meta.prepend(number);

        // move the size category (nibble / tidbit / chunk) out of the tag list
        post.querySelectorAll(".listing-category").forEach((cat) => {
            const name = cat.textContent.trim().toLowerCase();
            if (!(name in sizes)) return;
            const size = document.createElement("div");
            size.className = "listing-size";
            size.title = `size: ${name}`;
            size.innerHTML = [1, 2, 3]
                .map((k) => `<span class="pip${k <= sizes[name] ? " full" : ""}"></span>`)
                .join("") + `<span class="listing-size-name">${name}</span>`;
            meta.append(size);
            cat.style.display = "none";
        });

        const minutes = Number(post.dataset.listingReadingTimeSort);
        if (minutes > 0) {
            const reading = document.createElement("div");
            reading.className = "listing-reading";
            reading.textContent = `${minutes} min read`;
            meta.append(reading);
        }
    });
});

function activateTab(tab_title) {
  // sphinx-design tab buttons have role="tab" and match by their text content
  const inputs = document.querySelectorAll('.gui-tabs input');
  const labels = document.querySelectorAll('.gui-tabs label');

  // find the name of the input that we want to click
  for (const label of labels) {
    if (label.textContent.trim() == tab_title) {
      input_name = label.getAttribute('for');
      break;
    }
  }

  // find the input that needs to be clicked
  for (const inp of inputs) {
    if (inp.getAttribute('id').trim() == input_name) {
      inp.click();
      inp.parentElement.scrollIntoView({ behavior: 'smooth', block: 'start' });
      break;
    }
  }

}

document.addEventListener('DOMContentLoaded', () => {
  // Wire up image region clicks
  document.querySelectorAll('.gui-wrap a.region[data-tab]').forEach(el => {
    el.addEventListener('click', e => {
      e.preventDefault();
      activateTab(el.dataset.tab);
    });
  });
});
